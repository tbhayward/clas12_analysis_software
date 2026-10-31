#!/usr/bin/env python3
"""Reconstructed-only AAOgen -> DVCS leakage worker; no generated branches.
Run from clas12_analysis_software root. One HIPO per job, FIFO -> PhysicsEvents ROOT.
"""
import argparse, concurrent.futures, hashlib, json, os, pathlib, re, shlex, shutil, subprocess, threading
P=pathlib.Path
SCRIPTS={'pi0_as_dvcs':'processing_pi0_as_dvcs_stream.groovy'}
def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--config',default='processing_scripts/campaign_config.json');ap.add_argument('--node',type=int,required=True);ap.add_argument('--nodes',type=int,default=2);ap.add_argument('--workers',type=int,default=10);ap.add_argument('--periods',nargs='+',default=None);ap.add_argument('--only',nargs='+',choices=list(SCRIPTS),default=list(SCRIPTS));ap.add_argument('--max-events',type=int,default=10000);ap.add_argument('--output-root',default='rebuild_data_pilot_v1');ap.add_argument('--execute',action='store_true');ap.add_argument('--converter-timeout',type=int,default=180);ap.add_argument('--retry-stale',action='store_true',help='Remove stale FIFO/partial ROOT files only for jobs without a valid completion marker')
 a=ap.parse_args()
 if not 1<=a.node<=a.nodes or a.workers<1 or a.max_events<0:ap.error('invalid node, workers or max-events')
 root=P.cwd();cfg=json.loads(P(a.config).read_text());out=P(a.output_root);jobs=[];missing=[];seen=set()
 for j in cfg['jobs']:
  if j['mode']!='pi0_rec' or 'pi0_as_dvcs' not in a.only or (a.periods and j['period'] not in a.periods):continue
  inp=P(j['input'])
  if not inp.is_dir():missing.append(str(inp));continue
  files=sorted(inp.rglob('*.hipo'))
  if not files:missing.append(str(inp)+' (empty)');continue
  for f in files:
   key=('pi0_as_dvcs',str(f.resolve()))
   if key in seen:continue
   seen.add(key)
   h=hashlib.sha256(('pi0_as_dvcs'+'\0'+str(f.resolve())).encode()).digest()
   if int.from_bytes(h[:8],'big')%a.nodes+1!=a.node:continue
   dest=out/j['period']/'pi0_as_dvcs'/'parts'/(f.stem+'_'+hashlib.sha256(str(f.resolve()).encode()).hexdigest()[:10]+'.root')
   jobs.append((j,f,dest))
 print(f'[PLAN] node={a.node}/{a.nodes} files={len(jobs)} workers={a.workers} execute={a.execute}',flush=True)
 if missing:print('[MISSING_INPUTS]',*missing,sep='\n ',flush=True)
 if not a.execute:return
 if missing:raise SystemExit('Missing input directories: refusing partial run')
 if not jobs:raise SystemExit('No jobs on this node')
 for x in ('g++','root-config'):
  if not shutil.which(x):raise SystemExit('Missing '+x)
 jar=root/'processing_classes/dist/processing_classes.jar'
 if not jar.is_file():raise SystemExit('Missing '+str(jar))
 source=root/'processing_scripts/rebuild_pilot/convert_data_stream.cpp';converter=root/'processing_scripts/rebuild_pilot/convert_data_stream'
 validator_source=root/'processing_scripts/rebuild_pilot/validate_reconstructed_root.cpp';validator=root/'processing_scripts/rebuild_pilot/validate_reconstructed_root'
 schema_header=root/'processing_scripts/rebuild_pilot/reconstructed_data_schema.h'
 for required_source in (source,validator_source,schema_header):
  if not required_source.is_file():raise SystemExit('Missing '+str(required_source))
 flags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True));libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
 subprocess.run(['g++','-O3','-std=c++17',*flags,str(source),'-o',str(converter),*libs],check=True)
 subprocess.run(['g++','-O3','-std=c++17',*flags,str(validator_source),'-o',str(validator),*libs],check=True)
 counts={'done':0,'skipped':0};lock=threading.Lock()
 def run(job):
  j,hipo,dest=job;dest.parent.mkdir(parents=True,exist_ok=True)
  ok=P(str(dest)+'.ok')
  required={'input':str(hipo),'mode':'pi0_as_dvcs','period':j['period'],'max_events':a.max_events,'input_size':hipo.stat().st_size}
  fifo=P(str(dest)+'.events.tsv');tmp=P(str(dest)+'.partial.root');gl=P(str(dest)+'.processing.log');cl=P(str(dest)+'.converter.log')
  valid=False
  if dest.is_file() and dest.stat().st_size>0 and ok.is_file():
   try:
    marker=json.loads(ok.read_text())
    valid=all(marker.get(k)==v for k,v in required.items())
   except (ValueError,OSError):pass
  if valid:
   valid=subprocess.run([str(validator),str(dest),'pi0_as_dvcs',str(marker['selected_rows'])],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).returncode==0
  if valid:
   with lock:counts['done']+=1;counts['skipped']+=1
   print('[SKIP_VALID_MARKER]','pi0_as_dvcs',hipo.name,flush=True);return
  stale=[f for f in (dest,ok,fifo,tmp) if f.exists() or f.is_symlink()]
  if stale:
   if not a.retry_stale:raise RuntimeError('Stale FIFO or partial output (use --retry-stale after confirming no other data pilot is active): '+str(dest))
   # Only clean artifacts belonging to this incomplete job, never completed jobs.
   for f in stale:
    if f.is_dir():raise RuntimeError('Unexpected directory at output path: '+str(f))
    f.unlink()
   print('[RETRY_STALE]','pi0_as_dvcs',hipo.name,flush=True)
  os.mkfifo(fifo)
  try:
   with gl.open('w') as g,cl.open('w') as c:
    conv=subprocess.Popen([str(converter),str(fifo),str(tmp),'pi0_as_dvcs'],stdout=c,stderr=subprocess.STDOUT)
    script=root/'processing_scripts/rebuild_pilot'/SCRIPTS['pi0_as_dvcs']
    args=[str(root/'coatjava/bin/run-groovy'),'-cp',str(jar),str(script),str(hipo),str(fifo),'1',str(j['energy']),'11','1',str(a.max_events)]
    cmd='source source_file.txt; if ($status != 0) exit 90; module load qadb/3.4.1; if ($status != 0) exit 91; exec '+' '.join(shlex.quote(x) for x in args)
    try:
     proc=subprocess.run(['/bin/tcsh','-c',cmd],cwd=root,stdout=g,stderr=subprocess.STDOUT)
     if proc.returncode:conv.terminate()
     try:code=conv.wait(timeout=a.converter_timeout)
     except subprocess.TimeoutExpired:conv.kill();code=conv.wait()
    except BaseException:conv.kill();conv.wait();raise
   if proc.returncode or code or not tmp.is_file() or tmp.stat().st_size==0:raise RuntimeError(f'Processor={proc.returncode} converter={code}; see {gl} {cl}')
   log=gl.read_text(errors='replace');cv=cl.read_text(errors='replace')
   totals=re.findall(r'\[HIPO_FILE_DONE\].*?events_total=(\d+)',log);rows=re.findall(r'\[DATA_ROOT_ROWS\]\s+(\d+)',cv)
   if len(totals)!=1 or len(rows)!=1 or int(totals[0])<int(rows[0]) or (a.max_events>0 and int(totals[0])>a.max_events):raise RuntimeError(f'Invalid completion counts: processed={totals} selected={rows}; {hipo}')
   checked=subprocess.run([str(validator),str(tmp),'pi0_as_dvcs',rows[0]],capture_output=True,text=True)
   if checked.returncode:raise RuntimeError(f'ROOT schema validation failed for {hipo}: {checked.stderr.strip()}')
   # Data processors emit only selected candidates: selected rows <= processed events.
   os.replace(tmp,dest);ok.write_text(json.dumps({**required,'processed_events':int(totals[0]),'selected_rows':int(rows[0])})+'\n')
   with lock:counts['done']+=1
   print('[DONE]','pi0_as_dvcs',hipo.name,'events',totals[0],'selected',rows[0],flush=True)
  finally:
   if fifo.exists():fifo.unlink()
   if tmp.exists():tmp.unlink()
 failures=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  for f in concurrent.futures.as_completed([pool.submit(run,j) for j in jobs]):
   try:f.result()
   except Exception as exc:failures.append(str(exc));print('[FAILED]',exc,flush=True)
 print(f'[FINAL] node={a.node} completed={counts["done"]}/{len(jobs)} skipped={counts["skipped"]} failed={len(failures)}',flush=True)
 if failures:raise SystemExit(1)
if __name__=='__main__':main()
