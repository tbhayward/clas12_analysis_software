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
 a.add_argument('--config',default='processing_scripts/campaign_config.json');a.add_argument('--node',type=int,required=True);a.add_argument('--nodes',type=int,default=2);a.add_argument('--workers',type=int,default=10);a.add_argument('--only',nargs='+',choices=MODE);a.add_argument('--periods',nargs='+');a.add_argument('--output-root',default=None);a.add_argument('--execute',action='store_true');a.add_argument('--max-events',type=int,default=0);a.add_argument('--report-seconds',type=int,default=60)
 o=a.parse_args()
 if o.node<1 or o.node>o.nodes or o.workers<1: a.error('Invalid node/workers')
 root=P.cwd();scripts=root/'processing_scripts/rebuild_pilot';cfg=json.loads(P(o.config).read_text())
 out=P(o.output_root or (cfg['output_root']+'_paired_rebuild'))
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
 state={'done':0,'failed':0,'active':{}};started=time.monotonic();stop=threading.Event()
 def report():
  while not stop.wait(o.report_seconds):
   with lock: active=list(state['active'].items());done=state['done'];failed=state['failed']
   counts=[]
   for key,log in active:
    try:
     with open(log,'rb') as h:h.seek(max(0,os.fstat(h.fileno()).st_size-8192));tail=h.read().decode(errors='replace')
     matches=re.findall(r'\[HIPO_PROGRESS\].*?events_processed=(\d+)',tail)
     if matches:counts.append(f'{key}: {matches[-1]}')
    except OSError:pass
   print(f'[PROGRESS] node={o.node} complete={done}/{len(jobs)} active={len(active)} failed={failed} elapsed_min={(time.monotonic()-started)/60:.1f}'+(' | '+ ' | '.join(counts[:4]) if counts else ''),flush=True)
 def run(job):
  sample,hipo,destination,energy=job;destination.parent.mkdir(parents=True,exist_ok=True)
  if destination.exists():
   if destination.stat().st_size>0 and P(str(destination)+'.ok').is_file():
    with lock:state['done']+=1
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
     try:code=conv.wait(timeout=180)
     except subprocess.TimeoutExpired:
      print(f'[CONVERTER_TIMEOUT] {hipo}: converter did not exit within 180 seconds after Groovy. Killing converter.',flush=True)
      conv.kill();code=conv.wait()
    except BaseException:
     conv.kill();conv.wait();raise
   if proc.returncode or code or not tmp.is_file() or tmp.stat().st_size==0:raise RuntimeError(f'Failed {hipo}: Groovy={proc.returncode} converter={code}; logs {log} {clog}')
   os.replace(tmp,destination)
   P(str(destination)+'.ok').write_text(json.dumps({'input':str(hipo),'mode':sample[1],'period':sample[0],'max_events':o.max_events})+'\n')
   with lock:state['done']+=1
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
