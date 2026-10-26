#!/usr/bin/env python3
"""Two-node DVCS campaign: 5 required modes, nine workers/node, progress by input bytes.
Run --prepare-only once, then --execute --node 1/2 and --node 2/2 on separate hosts.
Byte-weighted progress is an estimate of workload, NOT an exact HIPO event percentage.
"""
import argparse,csv,json,os,re,shlex,subprocess,threading,time,shutil
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed

SCRIPTS={'dvcs_data':'processing_dvcs.groovy','pi0_data':'processing_exclusive_pi0.groovy','dvcs_rec':'processing_mc_rec_dvcs.groovy','pi0_rec':'processing_mc_rec_exclusive_pi0.groovy','pi0_as_dvcs':'processing_mc_rec_dvcs.groovy'}
lock=threading.Lock()

def prepare():
    jar=Path('processing_classes/dist/processing_classes.jar'); service=Path('coatjava/lib/services/processing_classes.jar')
    if not jar.is_file(): raise RuntimeError('Missing JAR; rebuild processing_classes first: '+str(jar))
    service.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(jar,service)
    for name in ('convert_txt_to_root','convert_dvcs_mc_diagnostics'):
        src=Path('processing_scripts')/(name+'.cpp')
        cflags=shlex.split(subprocess.check_output(['root-config','--cflags'],text=True))
        libs=shlex.split(subprocess.check_output(['root-config','--libs'],text=True))
        subprocess.run(['g++',*cflags,str(src),'-o',str(src.with_suffix('')),*libs],check=True)
    print('Preparation complete: JAR installed and converters compiled. No git pull performed.',flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',default='processing_scripts/campaign_config.json')
    ap.add_argument('--node',default='1/2',help='Node partition 1/2 or 2/2 (or 1/1)')
    ap.add_argument('--workers',type=int,default=9)
    ap.add_argument('--execute',action='store_true');ap.add_argument('--prepare-only',action='store_true')
    ap.add_argument('--skip-existing',action='store_true');ap.add_argument('--manifest',default=None)
    ap.add_argument('--output-root',default=None)
    ap.add_argument('--periods',nargs='*');ap.add_argument('--only',nargs='*',choices=sorted(SCRIPTS))
    ap.add_argument('--progress-seconds',type=int,default=60)
    a=ap.parse_args()
    if a.prepare_only: prepare();return
    if a.workers<1 or a.workers>9:ap.error('workers must be 1..9')
    m=re.fullmatch(r'(\d+)/(\d+)',a.node)
    if not m:ap.error('node must be N/T, e.g. 1/2')
    node,total=map(int,m.groups())
    if not 1<=node<=total:ap.error('invalid node partition')
    cfg=json.loads(Path(a.config).read_text()); root=Path(a.output_root or (cfg['output_root']+'_no_proton_eloss'))/'production'
    jobs=[]; missing=[]
    for j in cfg['jobs']:
        if j['mode'] not in SCRIPTS:continue
        if a.periods and j['period'] not in a.periods:continue
        if a.only and j['mode'] not in a.only:continue
        inp=Path(j['input'])
        if not inp.is_dir():missing.append(str(inp));continue
        files=sorted(inp.rglob('*.hipo'))
        if not files:missing.append(str(inp)+' (no HIPO files)');continue
        size=sum(f.stat().st_size for f in files)
        jobs.append(dict(j,files=len(files),bytes=size,source=str(inp.resolve())))
    if missing:
        print('Missing input directories:',*missing,sep='\n  ')
        if a.execute:raise SystemExit('Refusing incomplete production campaign')
    # Keep all modes reading a given directory on one node, and sequential within that node.
    groups=defaultdict(list)
    for j in jobs:groups[j['source']].append(j)
    bins=[{'size':0,'groups':[]} for _ in range(total)]
    for source,group in sorted(groups.items(),key=lambda kv:-sum(x['bytes'] for x in kv[1])):
        k=min(range(total),key=lambda i:bins[i]['size'])
        bins[k]['groups'].append((source,group));bins[k]['size']+=sum(x['bytes'] for x in group)
    selected=bins[node-1]['groups']; rows=[]
    for source,group in selected:
        for j in group:
            mode=j['mode']; current=j.get('current','nominal'); period=j['period']
            out=root/period/mode/current/'all'
            mc=mode not in ('dvcs_data','pi0_data')
            args=['csh','processing_scripts/processing.csh','processing_scripts/'+SCRIPTS[mode],j['input'],str(out),'0',str(j['energy']),'11' if mc else '0','1' if mc else '0','0']
            if mc:args.append(period)
            expected=[str(out)+'.root']+([str(out)+'.mc_events.root',str(out)+'.mc_particles.root'] if mc else [])
            rows.append(dict(node=a.node,period=period,mode=mode,current=current,input=j['input'],input_files=j['files'],input_bytes=j['bytes'],output=str(out),expected=';'.join(expected),command=shlex.join(args),source=source))
    manifest=Path(a.manifest or f'processing_scripts/campaign_manifest_node{node}of{total}.csv')
    manifest.parent.mkdir(parents=True,exist_ok=True)
    with manifest.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['node','period','mode','current','input','input_files','input_bytes','output','expected','command','source']);w.writeheader();w.writerows(rows)
    print(f'Node {a.node}: {len(rows)} jobs, {len(selected)} independent input directories, {sum(r["input_files"] for r in rows)} HIPO-file reads, {sum(r["input_bytes"] for r in rows)/1e12:.3f} TB scheduled input (counting repeated reads); workers={a.workers}; execute={a.execute}; manifest={manifest}',flush=True)
    for k,b in enumerate(bins,1):print(f'  node {k}/{total}: {sum(len(g) for _,g in b["groups"])} jobs; {b["size"]/1e12:.3f} TB input',flush=True)
    if not a.execute:return
    if not rows:raise SystemExit('No jobs selected')
    if not Path('processing_scripts/convert_txt_to_root').is_file() or not Path('processing_scripts/convert_dvcs_mc_diagnostics').is_file():raise SystemExit('Converters missing. Run --prepare-only once first.')
    for r in rows:
        for p in r['expected'].split(';'):
            if Path(p).exists() and not (a.skip_existing and all(Path(q).is_file() and Path(q).stat().st_size>0 for q in r['expected'].split(';'))):raise SystemExit('Existing output: '+p)
    state={'done':0,'bytes':0,'active':{},'fail':0};start=time.monotonic(); total_bytes=sum(r['input_bytes'] for r in rows)
    def run_group(group):
        for r in group:
            exp=[Path(p) for p in r['expected'].split(';')]
            if a.skip_existing and all(p.is_file() and p.stat().st_size>0 for p in exp):
                with lock:state['done']+=1;state['bytes']+=r['input_bytes']
                continue
            out=Path(r['output']);out.parent.mkdir(parents=True,exist_ok=True)
            log=Path(str(out)+'.processing.log')
            with lock:state['active'][r['output']]=dict(log=str(log),label=f'{r["period"]} {r["mode"]} {r["current"]}')
            env=os.environ.copy();env['DVCS_CAMPAIGN_WORKER']='1'
            with log.open('w') as f:res=subprocess.run(shlex.split(r['command']),stdout=f,stderr=subprocess.STDOUT,env=env)
            with lock:state['active'].pop(r['output'],None)
            if res.returncode or not all(p.is_file() and p.stat().st_size>0 for p in exp):
                with lock:state['fail']+=1
                raise RuntimeError(f'Failed {r["period"]} {r["mode"]} {r["current"]}, exit={res.returncode}; log={log}')
            with lock:state['done']+=1;state['bytes']+=r['input_bytes']
            print(f'[DONE] {r["period"]} {r["mode"]} {r["current"]} ({state["done"]}/{len(rows)} jobs)',flush=True)
    stop=threading.Event()
    def reporter():
        while not stop.wait(a.progress_seconds):
            with lock:done=state['done'];done_bytes=state['bytes'];active=list(state['active'].values());fail=state['fail']
            # The only defensible exact progress measure without pre-scanning HIPO is completed jobs/files.
            pct=100*done_bytes/total_bytes if total_bytes else 100
            samples=[]
            for item in active[:5]:
                try:
                    with open(item['log'],'rb') as f:
                        f.seek(max(0,os.fstat(f.fileno()).st_size-12000));tail=f.read().decode(errors='replace')
                    matches=re.findall(r'processed:\s*([\d,]+)\s*events',tail)
                    if matches:samples.append(item['label']+': '+matches[-1]+' events')
                except OSError:pass
            elapsed=(time.monotonic()-start)/60
            print(f'[PROGRESS] completed {done}/{len(rows)} jobs; {pct:.1f}% of scheduled input bytes in completed jobs; active={len(active)}; failures={fail}; elapsed={elapsed:.1f} min'+(' | '+ ' | '.join(samples) if samples else ''),flush=True)
    thread=threading.Thread(target=reporter,daemon=True);thread.start();errors=[]
    try:
        with ThreadPoolExecutor(max_workers=a.workers) as pool:
            futures={pool.submit(run_group,group):source for source,group in selected}
            for future in as_completed(futures):
                try:future.result()
                except Exception as e:errors.append(str(e));print('[ERROR]',e,flush=True)
    finally:stop.set();thread.join(timeout=2)
    if errors:raise SystemExit(f'{len(errors)} groups failed; inspect logs')
    print(f'COMPLETE node {a.node}: {len(rows)} jobs; elapsed={(time.monotonic()-start)/60:.1f} min',flush=True)
if __name__=='__main__':main()
