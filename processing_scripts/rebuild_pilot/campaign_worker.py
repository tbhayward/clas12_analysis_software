#!/usr/bin/env python3
"""Per-HIPO, two-node paired-MC campaign runner (FIFO path fix). Does not touch legacy outputs.
Supports dvcs_mc and pi0_mc; other modes intentionally rejected pending validation.
Run from clas12_analysis_software root. Both nodes use identical config and --nodes 2.
"""
import argparse, concurrent.futures, hashlib, json, os, pathlib, re, shlex, shutil, subprocess, threading, time
P=pathlib.Path
MODE={'dvcs_mc':'processing_mc_rec_dvcs.groovy','pi0_mc':'processing_mc_rec_exclusive_pi0.groovy'}
LEGACY={'dvcs_rec':'dvcs_mc','pi0_rec':'pi0_mc'}
lock=threading.Lock()
def main():
 a=argparse.ArgumentParser(description=__doc__)
 a.add_argument('--config',default='processing_scripts/campaign_config.json');a.add_argument('--node',type=int,required=True);a.add_argument('--nodes',type=int,default=2);a.add_argument('--workers',type=int,default=10);a.add_argument('--only',nargs='+',choices=MODE);a.add_argument('--periods',nargs='+');a.add_argument('--output-root',default=None);a.add_argument('--execute',action='store_true');a.add_argument('--max-events',type=int,default=0);a.add_argument('--report-seconds',type=int,default=60);a.add_argument('--converter-timeout',type=int,default=180)
 o=a.parse_args()
 if o.node<1 or o.node>o.nodes or o.workers<1: a.error('Invalid node/workers')
 root=P.cwd();scripts=root/'processing_scripts/rebuild_pilot';cfg=json.loads(P(o.config).read_text())
 out=P(o.output_root or (cfg['output_root']+'_paired_rebuild'))
 if o.max_events < 0:a.error('--max-events must be >= 0')
 if o.converter_timeout < 1:a.error('--converter-timeout must be positive')
 jobs=[];seen=set();missing=[]
 for j in cfg['jobs']:
  mode=LEGACY.get(j['mode']);
  if mode is None or (o.only and mode not in o.only) or (o.periods and j['period'] not in o.periods):continue
  inp=P(j['input'])
  if not inp.is_dir():missing.append(str(inp));continue
  files=sorted(inp.rglob('*.hipo'))
  if not files:missing.append(str(inp)+' (empty)');continue
  for f in files:
   key=(mode,str(f.resolve()))
   if key in seen:continue
   seen.add(key)
   digest=hashlib.sha256((mode+'\0'+str(f.resolve())).encode()).digest()
   node=int.from_bytes(digest[:8],'big')%o.nodes+1
   if node!=o.node:continue
   sample=(j['period'],mode,j.get('current','nominal'))
   destination=out.joinpath(*sample,'parts',f.stem+'_'+hashlib.sha256(str(f.resolve()).encode()).hexdigest()[:10]+'.root')
   jobs.append((sample,f,destination,j['energy']))
 print(f'[PLAN] node={o.node}/{o.nodes} files={len(jobs)} workers={o.workers} execute={o.execute}',flush=True)
 if missing:print('[MISSING INPUTS]',*missing,sep='\n  ',flush=True)
 if not o.execute:return
 if missing:raise SystemExit('Missing inputs: refusing incomplete campaign')
 if not jobs:raise SystemExit('No jobs for this node')
 for x in ('root-config','g++'):
  if not shutil.which(x):raise SystemExit('Missing '+x)
 jar=root/'processing_classes/dist/processing_classes.jar'
 if not jar.is_file():raise SystemExit('Missing '+str(jar))
 converter=scripts/'convert_paired_stream'
 flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True));libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
 subprocess.run(['g++','-O3','-std=c++17',*flags,str(scripts/'convert_paired_stream.cpp'),'-o',str(converter),*libs],check=True)
 state={'done':0,'failed':0,'active':{},'events_done':0};started=time.monotonic();stop=threading.Event()
 def report():
  while not stop.wait(o.report_seconds):
   with lock: active=list(state['active'].items());done=state['done'];failed=state['failed'];finished_events=state['events_done']
   counts=[]; processed=finished_events; expected=finished_events
   for key,log in active:
    try:
     with open(log,'rb') as h:h.seek(max(0,os.fstat(h.fileno()).st_size-8192));tail=h.read().decode(errors='replace')
     matches=re.findall(r'\[HIPO_PROGRESS\].*?events_processed=(\d+)',tail)
     if matches:counts.append(f'{key}: {matches[-1]}')
     processed+=int(matches[-1]) if matches else 0
     total=re.findall(r'number of\s+events\s*:\s*(\d+)',tail)
     if total:expected+=min(int(total[-1]),o.max_events) if o.max_events else int(total[-1])
    except OSError:pass
   elapsed=time.monotonic()-started
   rate=processed/elapsed if elapsed>0 else 0
   print(f'[PROGRESS] node={o.node} complete={done}/{len(jobs)} active={len(active)} failed={failed} events_observed={processed} rate_ev_s={rate:.1f} elapsed_min={elapsed/60:.1f}'+(' | '+ ' | '.join(counts[:4]) if counts else ''),flush=True)
 def run(job):
  sample,hipo,destination,energy=job;destination.parent.mkdir(parents=True,exist_ok=True)
  if destination.exists():
   if destination.stat().st_size>0 and P(str(destination)+'.ok').is_file():
    try: marker=json.loads(P(str(destination)+'.ok').read_text())
    except (OSError,ValueError) as exc:raise RuntimeError('Unreadable completion marker '+str(destination)+': '+str(exc))
    required={'input':str(hipo),'mode':sample[1],'period':sample[0],'max_events':o.max_events}
    if any(marker.get(k)!=v for k,v in required.items()):raise RuntimeError('Output belongs to another run/limit; use a different --output-root: '+str(destination))
    if marker.get('input_size') is not None and marker['input_size']!=hipo.stat().st_size:raise RuntimeError('Input HIPO size changed: '+str(hipo))
    with lock:state['done']+=1;state['events_done']+=int(marker.get('processed_events',0))
    print('[SKIP_VALID_MARKER]',sample,hipo.name,flush=True)
    return
   raise RuntimeError('Existing incomplete output, refusing overwrite: '+str(destination))
  selected=P(str(destination)+'.selected.txt');fifo=P(str(selected)+'.events.tsv');log=P(str(destination)+'.processing.log');clog=P(str(destination)+'.converter.log');tmp=P(str(destination)+'.partial.root')
  for p in (fifo,tmp):
   if p.exists():raise RuntimeError('Stale partial/FIFO: '+str(p))
  if P(str(selected)+'.events.tsv').exists():raise RuntimeError('Stream path already exists: '+str(selected)+'.events.tsv')
  os.mkfifo(fifo)
  with lock:state['active'][hipo.name]=str(log)
  try:
   with log.open('w') as gl,clog.open('w') as cl:
    conv=subprocess.Popen([str(converter),str(fifo),str(tmp)],stdout=cl,stderr=subprocess.STDOUT)
    args=[str(root/'coatjava/bin/run-groovy'),'-cp',str(jar),str(scripts/MODE[sample[1]]),str(hipo),str(selected),'1',str(energy),'11','1',str(o.max_events),sample[0]]
    command='source source_file.txt; if ($status != 0) exit 90; module load qadb/3.4.1; if ($status != 0) exit 91; exec '+' '.join(shlex.quote(s) for s in args)
    try:
     proc=subprocess.run(['/bin/tcsh','-c',command],cwd=root,stdout=gl,stderr=subprocess.STDOUT)
     if proc.returncode:conv.terminate()
     try:code=conv.wait(timeout=o.converter_timeout)
     except subprocess.TimeoutExpired:
      print(f'[CONVERTER_TIMEOUT] {hipo}: converter did not exit within {o.converter_timeout} seconds after Groovy. Killing converter.',flush=True)
      conv.kill();code=conv.wait()
    except BaseException:
     conv.kill();conv.wait();raise
   if proc.returncode or code or not tmp.is_file() or tmp.stat().st_size==0:raise RuntimeError(f'Failed {hipo}: Groovy={proc.returncode} converter={code}; logs {log} {clog}')
   # Check converter row count against Groovy completion count before accepting output.
   converter_log=clog.read_text(errors='replace')
   groovy_log=log.read_text(errors='replace')
   rows=re.findall(r'\[PAIRED_ROOT_ROWS\]\s+(\d+)',converter_log)
   totals=re.findall(r'\[HIPO_FILE_DONE\].*?events_total=(\d+)',groovy_log)
   # Some Groovy versions emit HIPO_FILE_DONE for each event rather than once.
   # Accept strictly increasing progress markers only if the last count matches
   # the converter's single final row count.
   counts=[int(value) for value in totals]
   increasing=all(previous < current for previous,current in zip(counts,counts[1:]))
   if len(rows)!=1 or not counts or not increasing or int(rows[0])!=counts[-1]:
    raise RuntimeError(f'Event-count mismatch or missing completion markers for {hipo}: converter={rows} Groovy={totals[-3:]}')
   if len(counts)>1:
    print(f'[LEGACY_PROGRESS_MARKERS] {hipo.name}: {len(counts)} progress markers; final count agrees with converter',flush=True)
   n=int(rows[0])
   if o.max_events and n>o.max_events:raise RuntimeError(f'Exceeded event cap for {hipo}: {n}>{o.max_events}')
   if n<1:raise RuntimeError('Empty ROOT event stream: '+str(hipo))
   os.replace(tmp,destination)
   marker={'input':str(hipo),'mode':sample[1],'period':sample[0],'max_events':o.max_events,'processed_events':n,'input_size':hipo.stat().st_size}
   P(str(destination)+'.ok').write_text(json.dumps(marker)+'\n')
   with lock:state['done']+=1;state['events_done']+=n
   print('[DONE]',sample,hipo.name,flush=True)
  finally:
   with lock:state['active'].pop(hipo.name,None)
   if fifo.exists():fifo.unlink()
 stopthread=threading.Thread(target=report,daemon=True);stopthread.start();fail=[]
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=o.workers) as pool:
   fs={pool.submit(run,j):j for j in jobs}
   for future in concurrent.futures.as_completed(fs):
    try:future.result()
    except Exception as exc:
     with lock:state['failed']+=1
     fail.append(str(exc));print('[FAILED]',exc,flush=True)
 finally:stop.set();stopthread.join(timeout=2)
 print(f'[FINAL] node={o.node} completed={state["done"]}/{len(jobs)} failed={len(fail)}',flush=True)
 if fail:raise SystemExit(1)
if __name__=='__main__':main()
