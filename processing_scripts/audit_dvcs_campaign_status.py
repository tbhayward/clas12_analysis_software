#!/usr/bin/env python3
"""Read-only audit of a running DVCS paired-MC campaign.
Run from clas12_analysis_software repository root. Python standard library only.
"""
import argparse, csv, os, re, sys
from pathlib import Path
from collections import defaultdict
from datetime import datetime

MODES = ('dvcs_data','pi0_data','dvcs_rec','pi0_rec','pi0_as_dvcs')
PERIODS = ('sp18_inb','sp18_out','fa18_inb','fa18_out','sp19_inb')

def fmtbytes(n):
    for unit in ('B','KB','MB','GB','TB'):
        if abs(n)<1000: return f'{n:.1f} {unit}'
        n/=1000
    return f'{n:.1f} PB'

def loginfo(path):
    if not path.is_file(): return None, '', ''
    try:
        with path.open('rb') as f:
            f.seek(max(0,path.stat().st_size-65536))
            tail=f.read().decode('utf-8','replace')
        matches=re.findall(r'processed:\s*([\d,]+)\s*events',tail,re.I)
        last=int(matches[-1].replace(',','')) if matches else None
        errors=re.findall(r'(?im)^.*(?:exception|traceback|outofmemoryerror|segmentation fault|error:).{0,150}$',tail)
        return last, datetime.fromtimestamp(path.stat().st_mtime).strftime('%m-%d %H:%M'), errors[-1].strip() if errors else ''
    except OSError as e: return None,'',str(e)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest-dir',default='processing_scripts')
    p.add_argument('--manifest',action='append',help='Explicit manifest CSV; may be repeated')
    p.add_argument('--output',default='processing_scripts/campaign_status_audit.csv')
    p.add_argument('--show-all',action='store_true',help='Show every sample, including completed ones')
    a=p.parse_args()
    paths=[Path(x) for x in a.manifest] if a.manifest else sorted(Path(a.manifest_dir).glob('campaign_manifest_node*of*.csv'))
    if not paths:
        print('ERROR: no campaign_manifest_node*of*.csv files found. Pass --manifest path/to/file.csv (repeat for both nodes).',file=sys.stderr)
        return 2
    print('MANIFESTS:', ', '.join(map(str,paths)))
    rows=[]; keys=set()
    for path in paths:
        with path.open(newline='') as f:
            for r in csv.DictReader(f):
                key=(r.get('period'),r.get('mode'),r.get('current'),r.get('output'))
                if key in keys: continue
                keys.add(key)
                output=Path(r['output']); expected=[Path(x) for x in r.get('expected','').split(';') if x]
                if not expected:
                    expected=[Path(str(output)+'.root')]
                    if r['mode'] not in ('dvcs_data','pi0_data'):
                        expected += [Path(str(output)+'.mc_events.root'),Path(str(output)+'.mc_particles.root')]
                sizes=[q.stat().st_size if q.is_file() else 0 for q in expected]
                ready=all(s>0 for s in sizes)
                log=Path(str(output)+'.processing.log')
                ev,stamp,error=loginfo(log)
                # Presence of nonzero ROOT files is NOT a reliable completed marker during conversion.
                # Treat it as outputs-present, not launcher-confirmed completion.
                txt=Path(str(output)+'.txt')
                tsv1=Path(str(output)+'.txt.events.tsv')
                tsv2=Path(str(output)+'.txt.particles.tsv')
                present=sum(q.is_file() for q in expected)
                if ready: status='OUTPUTS_PRESENT_VERIFY_DONE'
                elif present or txt.is_file() or tsv1.is_file() or tsv2.is_file(): status='IN_PROGRESS_OR_INTERRUPTED'
                elif log.is_file(): status='LOG_ONLY'
                else: status='NOT_STARTED'
                rows.append(dict(node=r.get('node','?'),period=r['period'],mode=r['mode'],current=r.get('current','nominal'),status=status,
                    input_bytes=int(r.get('input_bytes') or 0),input_files=int(r.get('input_files') or 0),
                    expected_count=len(expected),root_present=present,root_nonzero=sum(s>0 for s in sizes),root_bytes=sum(sizes),
                    expected_files=';'.join(str(q) for q in expected),missing_files=';'.join(str(q) for q,s in zip(expected,sizes) if s==0),
                    processed_events=ev if ev is not None else '',log_updated=stamp,log_error=error,
                    txt_bytes=txt.stat().st_size if txt.is_file() else 0,
                    events_tsv_bytes=tsv1.stat().st_size if tsv1.is_file() else 0,
                    particles_tsv_bytes=tsv2.stat().st_size if tsv2.is_file() else 0,
                    log_path=str(log),output=str(output)))
    rows.sort(key=lambda r:(PERIODS.index(r['period']) if r['period'] in PERIODS else 99,
                            MODES.index(r['mode']) if r['mode'] in MODES else 99,r['current']))
    print('\nSTATUS BY RUN PERIOD AND MODE (root file presence is provisional until launcher prints [DONE]):')
    for period in sorted(set(r['period'] for r in rows),key=lambda x:PERIODS.index(x) if x in PERIODS else 99):
        print('\n'+period.upper())
        for r in rows:
            if r['period']!=period:continue
            print(f"  {r['mode']:<12} {r['current']:<9} {r['node']:<4} {r['status']:<29} ROOT={r['root_nonzero']}/{r['expected_count']} ({fmtbytes(r['root_bytes']):>10}) input={fmtbytes(r['input_bytes']):>10} events={str(r['processed_events']):>11} log={r['log_updated']}")
            if r['log_error']:print('     LOG WARNING:',r['log_error'][:200])
    print('\nSUMMARY BY MODE:')
    for mode in MODES:
        subset=[r for r in rows if r['mode']==mode]
        if not subset:continue
        ready=[r for r in subset if r['status']=='OUTPUTS_PRESENT_VERIFY_DONE']
        total=sum(r['input_bytes'] for r in subset)
        print(f"  {mode:<12} outputs present {len(ready):>2}/{len(subset):<2}  input bytes associated with outputs present {100*sum(r['input_bytes'] for r in ready)/total if total else 0:5.1f}%  scheduled={fmtbytes(total)}")
    print('\nSUMMARY BY NODE:')
    for node in sorted(set(r['node'] for r in rows)):
        subset=[r for r in rows if r['node']==node]
        ready=[r for r in subset if r['status']=='OUTPUTS_PRESENT_VERIFY_DONE']
        total=sum(r['input_bytes'] for r in subset)
        print(f"  {node}: {len(ready)}/{len(subset)} outputs present, {100*sum(r['input_bytes'] for r in ready)/total if total else 0:.1f}% of scheduled bytes associated with these jobs; scheduled {fmtbytes(total)}")
    print('\nLARGEST JOBS WITHOUT ALL OUTPUTS PRESENT:')
    for r in sorted((r for r in rows if r['status']!='OUTPUTS_PRESENT_VERIFY_DONE'),key=lambda x:-x['input_bytes'])[:12]:
        print(f"  {r['period']} {r['mode']} {r['current']}: input={fmtbytes(r['input_bytes'])}, events={r['processed_events'] or '?'}, roots={r['root_nonzero']}/{r['expected_count']}, status={r['status']}")
    dest=Path(a.output);dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(f'\nWrote {len(rows)} rows to {dest}')
    print('IMPORTANT: launcher [DONE] messages, not file presence, establish successful completion. The launcher currently credits zero partial input bytes for running jobs.')
    return 0
if __name__=='__main__':sys.exit(main())
