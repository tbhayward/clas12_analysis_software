#!/usr/bin/env python3
"""Read-only report for an existing campaign; no workers are started or stopped."""
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import time

def read(path):
    try:return json.loads(path.read_text())
    except FileNotFoundError:return {}

def tail(path):
    try:
        with path.open('rb') as stream:
            stream.seek(0,2);size=stream.tell();stream.seek(max(0,size-65536))
            return size,stream.read().decode('utf-8',errors='replace')
    except FileNotFoundError:return 0,''

def report(root):
    manifest=read(root/'manifest.json')
    if not manifest:raise RuntimeError('No manifest found: '+str(root))
    print('Campaign:',root,'  observed on:',socket.gethostname())
    print('Batch                            files   GiB   elapsed_s process_s convert_s    events events/s log_MiB state')
    for job in manifest['jobs']:
        jid=job['id'];done=read(root/'chunks'/(jid+'.json'));owner=read(root/'claims'/jid/'owner.json');failed=read(root/'claims'/jid/'failed.json')
        log_bytes,text=tail(root/'logs'/(jid+'.log'))
        counters=done.get('diagnostics',{}).get('counts',{})
        if not counters:
            summaries=re.findall(r'SUMMARY[^\n]*',text)
            if summaries:counters={k:int(v) for k,v in re.findall(r'(\w+)=(\d+)',summaries[-1])}
        events=counters.get('seen');elapsed=done.get('seconds')
        if elapsed is None and owner:elapsed=max(0,time.time()-owner['started'])
        processing=done.get('processing_seconds');conversion=done.get('conversion_seconds')
        rate=events/processing if events is not None and processing and processing>0 else None
        fmt=lambda x:f'{x:.1f}' if x is not None else '-'
        state='DONE' if done else 'FAILED' if failed else 'ACTIVE' if owner else 'QUEUED'
        print(f'{jid:45} {len(job["files"]):3} {job["bytes"]/1024**3:6.3f} {fmt(elapsed):>9} {fmt(processing):>9} {fmt(conversion):>9} {str(events) if events is not None else "-":>9} {fmt(rate):>8} {log_bytes/1024**2:7.2f} {state}',flush=True)
        if log_bytes>10*1024**2:
            print('  Large log; last lines:',repr('\n'.join(text.splitlines()[-3:])))
    print('\nLocal Java/Groovy processes: PID PPID %CPU %MEM ELAPSED CPU_TIME COMMAND')
    result=subprocess.run(['ps','-u',str(os.getuid()),'-o','pid,ppid,pcpu,pmem,etime,time,args'],capture_output=True,text=True)
    if result.returncode==0:
        for line in result.stdout.splitlines()[1:]:
            if re.search(r'java|groovy',line,re.I):print(line)
    else:print(result.stderr.strip())
    processor=root/'build/source/processing_exclusive_piminus.groovy'
    if processor.exists():
        source=processor.read_text()
        print('\nFrozen processor QA before fitter:',source.find('qa.pass(')<source.find('fitter.select('))
    print('Read-only report; no campaign files or processes changed.')

if __name__=='__main__':
    if len(sys.argv)!=2:raise SystemExit('Usage: performance_report.py /path/to/campaign')
    report(Path(sys.argv[1]).resolve())
