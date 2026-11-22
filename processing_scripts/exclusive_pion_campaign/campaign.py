#!/usr/bin/env python3
"""Exclusive pion campaign. Start workers separately on each farm host; no SSH is used."""
import argparse, concurrent.futures, hashlib, json, os, re, shlex, shutil, socket, subprocess, sys, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
HOSTS=('ifarm2401','ifarm2402')
EXCLUDED={5247,5345,17331,16987,17079,17190,17639,16850,16851,16852,16855,16879}
def allowed(run): return run not in EXCLUDED and not 16600<run<16700 and not 17768<run<=17811

def datasets():
    out=[]
    for period,torus,version in [('spring2019','-1','v0'),('fall2019','+1','v1'),('spring2020','-1','v1')]:
        out.append(dict(id=f'RGB_{period}_D2',target='D2',path=f'/cache/clas12/rg-b/production/recon/{period}/torus{torus}/pass2/{version}/dst/train/sidisdvcs'))
    for period,short in [('summer2022','su22'),('fall2022','fa22'),('spring2023','sp23')]:
        base={'summer2022':'summer22/pass1/10.5gev','fall2022':'fall22/pass1','spring2023':'spring23/pass1'}[period]
        for target in ['ND3','C','CH2']+(['CD2'] if short=='sp23' else []):
            out.append(dict(id=f'RGC_{period}_{target}',target=target,path=f'/cache/clas12/rg-c/production/{base}/{target}/dst/train/sidisdvcs'))
        for target in ['ET','He']:
            out.append(dict(id=f'RGC_{period}_{target}',target=target,path=f'/work/clas12/thayward/CLAS12_SIDIS/rgc_ET_HE_runs/{short}/{target}'))
    return out

def atomic_json(path,value):
    tmp=path.with_name(path.name+f'.{os.getpid()}.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');os.replace(tmp,path)
def read(path): return json.loads(path.read_text())
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def run(cmd,**kwargs): return subprocess.run([str(x) for x in cmd],check=True,**kwargs)
def classnamepath(build):
    extra=read(build/'extra_classpath.json') if (build/'extra_classpath.json').exists() else []
    return os.pathsep.join([str(build/'classes')]+[str(p) for p in sorted((build/'lib').glob('*.jar'))]+[str(build/p) for p in extra]+[str(build/'qadb'),str(build/'qadb/src'),str(build/'qadb/src/main/groovy')])
def run_number(path):
    # Standard CLAS12 train names: ..._016123.hipo / ..._016123_00000.hipo.
    hits=re.findall(r'(?<!\d)(0*[1-9]\d{3,4})(?!\d)',path.stem)
    valid={int(x) for x in hits if 1000<=int(x)<=29999}
    if len(valid)==1:return valid.pop()
    # Explicit run-number parent directories are also supported.
    parents={int(p) for p in path.parts[:-1] if re.fullmatch(r'0*[1-9]\d{3,4}',p) and 1000<=int(p)<=29999}
    if len(parents)==1:return parents.pop()
    raise ValueError(f'Cannot unambiguously identify run for {path}; supply --run-map JSON mapping absolute paths to runs')
def make_batches(files,target_bytes):
    batches=[];current=[];size=0
    for f in sorted(files,key=lambda x:x['path']):
        if current and size+f['size']>target_bytes:batches.append(current);current=[];size=0
        current.append(f);size+=f['size']
    if current:batches.append(current)
    return batches

def prepare(a):
    root=a.campaign.resolve()
    root.mkdir(parents=True,exist_ok=True)
    lock=root/'prepare.claim';lock.mkdir()  # exclusive preparation; never overwrite a campaign
    try:
        if (root/'manifest.json').exists():raise RuntimeError('Campaign already prepared; use status/worker, or choose a new directory')
        if any(root.glob('build*')):raise RuntimeError('Partial preparation exists; choose a fresh campaign directory')
        repo=a.repo.resolve();inventory=[];jobs=[];overrides=read(a.run_map) if a.run_map else {};select=read(a.test_runs) if a.test_runs else {}
        for ds in datasets():
            source=Path(ds['path']);groups={}
            if not source.is_dir():raise RuntimeError(f'Missing input directory: {source}')
            for p in sorted(source.rglob('*.hipo')):
                stat=p.stat()
                if stat.st_size<=0:raise RuntimeError(f'Empty HIPO input: {p}')
                if not os.access(p,os.R_OK):raise RuntimeError(f'Unreadable input: {p}')
                n=int(overrides[str(p)]) if str(p) in overrides else run_number(p)
                groups.setdefault(n,[]).append(dict(path=str(p),size=stat.st_size,mtime_ns=stat.st_mtime_ns))
            eligible=sorted(n for n in groups if allowed(n))
            if not eligible:raise RuntimeError(f'No eligible input runs: {ds["id"]}')
            selected=eligible
            if a.mode=='test':
                # Prefer a small run for a quick smoke test; user can select alternatives.
                chosen=int(select[ds['id']]) if ds['id'] in select else min(eligible,key=lambda n:(sum(f['size'] for f in groups[n]),n))
                if chosen not in eligible:raise RuntimeError(f'Invalid test run {chosen} for {ds["id"]}')
                selected=[chosen]
            record=dict(ds,selected_runs=selected,all_runs=sorted(groups),excluded_runs=sorted(set(groups)-set(eligible)),files=groups)
            inventory.append(record)
            for n in selected:
                for ix,batch in enumerate(make_batches(groups[n],int(a.batch_gib*1024**3))):
                    jid=f'{ds["id"]}_r{n:06d}_b{ix:04d}'
                    jobs.append(dict(id=jid,dataset=ds['id'],target=ds['target'],run=n,files=batch,bytes=sum(f['size'] for f in batch)))
        # Freeze campaign sources, dependency jars, calibration helpers, CSV, and QA database.
        build=root/'build';shutil.copytree(HERE,build/'source',ignore=shutil.ignore_patterns('__pycache__','build'))
        (build/'lib').mkdir();(build/'classes').mkdir()
        jars=[repo/'processing_classes/dist/processing_classes.jar']
        for sub in ['clas','services','utils']:
            jars.extend(p for p in (repo/'coatjava/lib'/sub).glob('*.jar') if p.name!='processing_classes.jar')
        for p in jars:
            dest=build/'lib'/p.name
            if dest.exists() and digest(dest)!=digest(p):raise RuntimeError(f'Conflicting JAR names: {p}')
            shutil.copy2(p,dest)
        # Use the installed QADB module, ahead of any stale bundled Groovy source.
        qadir=os.environ.get('QADB')
        if not qadir or not Path(qadir).is_dir():raise RuntimeError('Load qadb/3.4.1 first (QADB must point to its database)')
        shutil.copytree(qadir,build/'qadb',symlinks=False)
        cp_entries=os.environ.get('CLASSPATH','').split(os.pathsep)
        extras=[]
        for entry in cp_entries:
            if not entry or entry=='.':continue
            p=Path(entry)
            if p.is_file() and p.suffix=='.jar':
                dest=build/'lib'/p.name
                if dest.exists() and digest(dest)!=digest(p):raise RuntimeError(f'Conflicting module JAR: {p}')
                shutil.copy2(p,dest)
            elif p.is_dir():
                try:relative=p.resolve().relative_to(Path(qadir).resolve());extras.append(str(Path('qadb')/relative))
                except ValueError:
                    dest=build/('module_cp_'+str(len(extras)));shutil.copytree(p,dest);extras.append(dest.name)
            elif entry.endswith('/*'):
                for jar in p.parent.glob('*.jar'):shutil.copy2(jar,build/'lib'/jar.name)
        atomic_json(build/'extra_classpath.json',extras)
        # Current helper sources override any stale classes in the original JAR.
        helper=build/'helpers';helper.mkdir()
        sources=[]
        for rel in ['extended_kinematic_fitters/analysis_fitter.java','extended_kinematic_fitters/generic_tests.java','extended_kinematic_fitters/fiducial_cuts.java','extended_kinematic_fitters/pid_cuts.java','extended_kinematic_fitters/energy_loss_corrections.java','analyzers/BeamEnergy.java']:
            src=repo/'processing_classes/src'/rel;dest=helper/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest);sources.append(dest)
        sources+=list((build/'source/src').glob('*.java'))
        run(['javac','--release','11','-cp',classnamepath(build),'-d',build/'classes',*sources])
        flags=shlex.split(subprocess.check_output(['root-config','--cflags','--libs'],text=True))
        run(['g++','-O3','-std=c++17',build/'source/src/convert_exclusive_piminus_txt_to_root.cpp','-o',build/'converter',*flags])
        shutil.copy2(repo/'analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv',build/'run_info.csv')
        converter_env=os.environ.copy();converter_env['PION_CONVERTER']=str(build/'converter')
        run([sys.executable,'-m','unittest','discover','-s',build/'source/tests','-p','test_converter.py','-v'],env=converter_env)
        # Resolve/compile QADB now, rather than discovering a missing module in 20 workers.
        env=os.environ.copy();env['QADB']=str(build/'qadb')
        run(['groovy','-cp',classnamepath(build),'-e',"def q=new clasqa.QADB('latest'); q.checkForDefect('Misc'); println('QADB ready')"],env=env)
        run(['groovy','-cp',classnamepath(build),'-e',"new GroovyClassLoader().parseClass(new File(args[0])); println('Campaign processor compiled')",str(build/'source/processing_exclusive_piminus.groovy')],env=env)
        for name in ['claims','chunks','work','logs','reports','root','hosts']: (root/name).mkdir()
        atomic_json(root/'inventory.json',inventory)
        atomic_json(root/'manifest.json',dict(mode=a.mode,workers_per_host=10,hosts=HOSTS,jobs=sorted(jobs,key=lambda j:(-j['bytes'],j['id'])),datasets=datasets(),created=time.time(),repo=str(repo)))
        frozen={str(p.relative_to(build)):digest(p) for p in build.rglob('*') if p.is_file()}
        atomic_json(root/'checksums.json',frozen)
        print(f'Prepared {a.mode}: {len(jobs)} batches, 19 datasets in {root}',flush=True)
        for ds in inventory:print(f'  {ds["id"]}: runs {ds["selected_runs"]}')
    finally:lock.rmdir()

def verified(root):
    build=root/'build'
    for name,sha in read(root/'checksums.json').items():
        if digest(build/name)!=sha:raise RuntimeError(f'Frozen software modified: {name}')
    return build

def progress_snapshot(root,manifest=None):
    manifest=manifest or read(root/'manifest.json');done=failed=active=0;entries=completed_bytes=0;details=[]
    total_bytes=sum(j['bytes'] for j in manifest['jobs'])
    for job in manifest['jobs']:
        jid=job['id'];completion=root/'chunks'/(jid+'.json');claim=root/'claims'/jid
        if completion.exists():
            info=read(completion);done+=1;entries+=info['entries'];completed_bytes+=job['bytes'];continue
        if (claim/'failed.json').exists():failed+=1;continue
        if not claim.exists():continue
        active+=1
        try:owner=read(claim/'owner.json')
        except FileNotFoundError:owner={}
        line=f'{jid} @ {owner.get("host","claiming")}'
        progress=root/'work'/jid/'progress.json'
        try:
            info=read(progress);counts=info['counts'];age=max(0,time.time()-info['updated'])
            line+=f': {info["stage"]}, {counts["seen"]:,} events, {counts["rows"]:,} rows, QA rejected {counts.get("qa_reject",0):,}, zero helicity {counts.get("zero_helicity",0):,}; update {age:.0f}s ago'
        except FileNotFoundError:line+=f': starting/processing ({max(0,time.time()-owner.get("started",time.time())):.0f}s elapsed)'
        details.append(line)
    total=len(manifest['jobs']);queued=total-done-failed-active
    print(f'[{time.strftime("%H:%M:%S")}] BOTH HOSTS: {done}/{total} batches done; {active} active; {queued} queued; {failed} failed; {entries:,} ROOT entries; {completed_bytes/1024**3:.2f}/{total_bytes/1024**3:.2f} GiB input completed',flush=True)
    for detail in details:print('  '+detail,flush=True)
    return dict(done=done,active=active,queued=queued,failed=failed,entries=entries)

def monitor(a):
    root=a.campaign.resolve();manifest=read(root/'manifest.json')
    while True:
        state=progress_snapshot(root,manifest)
        if a.once or state['done']+state['failed']==len(manifest['jobs']):break
        time.sleep(a.interval)

def process_job(root,build,j):
    claim=root/'claims'/j['id']
    try:claim.mkdir()
    except FileExistsError:return
    atomic_json(claim/'owner.json',dict(host=socket.gethostname(),pid=os.getpid(),started=time.time()))
    done=root/'chunks'/(j['id']+'.json')
    try:
        if done.exists():
            (claim/'owner.json').unlink();claim.rmdir();return
        work=root/'work'/j['id'];work.mkdir(exist_ok=True)
        for f in j['files']:
            s=Path(f['path']).stat()
            if s.st_size!=f['size'] or s.st_mtime_ns!=f['mtime_ns']:raise RuntimeError('Input changed: '+f['path'])
        listing=work/'inputs.list';listing.write_text(''.join(f['path']+'\n' for f in j['files']))
        text=work/'events.txt';tmp=work/'events.root';env=os.environ.copy();env['QADB']=str(build/'qadb');env['JAVA_OPTS']=os.environ.get('PION_JAVA_OPTS','-Dsun.java2d.pmoffscreen=false -Xms1024m -Xmx2048m')
        started=time.time()
        print(f'START {j["id"]} @ {socket.gethostname()}: {len(j["files"])} files, {j["bytes"]/1024**3:.2f} GiB',flush=True)
        with (root/'logs'/(j['id']+'.log')).open('w') as log:
            processor=build/'source/processing_exclusive_piminus.groovy'
            command=['groovy','-cp',classnamepath(build),processor,listing,text,str(j['run'])]
            # Earlier frozen processors accept exactly three arguments.
            if not processor.exists() or 'args.length!=3' not in processor.read_text():command.append(work/'progress.json')
            runtime=dict(java_opts=env['JAVA_OPTS'],runner_sha256=digest(Path(__file__).resolve()),started=started)
            atomic_json(work/'runtime.json',runtime)
            print(f'  Java settings: {env["JAVA_OPTS"]}',flush=True)
            run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
            processing_seconds=time.time()-started
            progress_path=work/'progress.json'
            if progress_path.exists():
                info=read(progress_path);info.update(stage='converting',updated=time.time());atomic_json(progress_path,info)
            print(f'CONVERT {j["id"]}: processing finished in {processing_seconds:.1f}s',flush=True)
            result=subprocess.run([str(x) for x in [build/'converter',text,tmp,build/'run_info.csv',j['target'],str(j['run'])]],capture_output=True,text=True)
            log.write(result.stderr)
            if result.returncode:raise RuntimeError('Converter: '+result.stderr.strip())
            entries=int(result.stdout.strip().splitlines()[-1])
            final=root/'chunks'/(j['id']+'.root');os.replace(tmp,final)
            diagnostic=read(work/'progress.json') if (work/'progress.json').exists() else {}
            if diagnostic:
                diagnostic.update(stage='done',updated=time.time());atomic_json(work/'progress.json',diagnostic)
            atomic_json(done,dict(runtime=runtime,diagnostics=diagnostic,entries=entries,seconds=time.time()-started,processing_seconds=processing_seconds,conversion_seconds=time.time()-started-processing_seconds,host=socket.gethostname(),bytes=j['bytes'],root_sha256=digest(final)))
        text.unlink();(claim/'failed.json').unlink(missing_ok=True)
        print(f'DONE {j["id"]}: {entries:,} entries; processing {processing_seconds:.1f}s, conversion {time.time()-started-processing_seconds:.1f}s',flush=True)
        if entries==0:
            print('  EMPTY TREE diagnostics: '+json.dumps(diagnostic.get('counts',{}),sort_keys=True),flush=True)
    except Exception as e:
        atomic_json(claim/'failed.json',dict(error=str(e),time=time.time()))
        print(f'FAILED {j["id"]}: {e}',flush=True)
        # Keep failure claim until an explicit retry, so another host cannot loop on it.
        return
    (claim/'owner.json').unlink(missing_ok=True);claim.rmdir()

def worker(a):
    root=a.campaign.resolve()
    if (root/'claims'/'conversion_recovery').exists():raise RuntimeError('Conversion recovery is running; do not start workers')
    manifest=read(root/'manifest.json');host=socket.gethostname().split('.')[0]
    if host not in HOSTS and not a.allow_other_host:raise RuntimeError('Workers are restricted to ifarm2401/ifarm2402')
    hostlock=root/'hosts'/host;hostlock.mkdir()
    try:
        atomic_json(hostlock/'owner.json',dict(pid=os.getpid(),started=time.time()))
        build=verified(root)
        print(f'{host}: starting up to {a.workers} local workers',flush=True)
        last_report=0.0
        with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
            # Refill from the shared queue. A claimed batch is never considered available.
            pending=set()
            while True:
                available=[j for j in manifest['jobs'] if not (root/'chunks'/(j['id']+'.json')).exists() and not (root/'claims'/j['id']).exists()]
                while available and len(pending)<a.workers:pending.add(pool.submit(process_job,root,build,available.pop(0)))
                if not pending:break
                if time.monotonic()-last_report>=30:
                    progress_snapshot(root,manifest);last_report=time.monotonic()
                finished,pending=concurrent.futures.wait(pending,timeout=1,return_when=concurrent.futures.FIRST_COMPLETED)
                for f in finished:f.result()
        complete=status(root)
        if complete:
            finalizer=root/'claims'/'finalize'
            try:finalizer.mkdir()
            except FileExistsError:print('Another host is merging the completed campaign',flush=True)
            else:
                try:merge(a)
                finally:finalizer.rmdir()
        else:
            failures=[j for j in manifest['jobs'] if (root/'claims'/j['id']/'failed.json').exists()]
            if failures:
                for job in failures[:3]:
                    failure=read(root/'claims'/job['id']/'failed.json')
                    print(f'FAILED {job["id"]}: {failure["error"]}; log: {root/"logs"/(job["id"]+".log")}',flush=True)
                raise RuntimeError(f'{len(failures)} failed batches; inspect logs before retrying')
            print('Local workers finished; remaining batches are claimed by other workers.',flush=True)
    finally:
        (hostlock/'owner.json').unlink(missing_ok=True);hostlock.rmdir()

def status(root):
    m=read(root/'manifest.json');total=0;failed=0
    for ds in m['datasets']:
        jobs=[j for j in m['jobs'] if j['dataset']==ds['id']];n=sum((root/'chunks'/(j['id']+'.json')).exists() for j in jobs);f=sum((root/'claims'/j['id']/'failed.json').exists() for j in jobs);total+=n;failed+=f
        print(f'{ds["id"]}: {n}/{len(jobs)} batches, {f} failed',flush=True)
    print(f'Total: {total}/{len(m["jobs"])} complete; {failed} failed',flush=True)
    return total==len(m['jobs'])

def merge(a):
    root=a.campaign.resolve();build=verified(root);m=read(root/'manifest.json');fail=False
    for ds in m['datasets']:
        jobs=sorted((j for j in m['jobs'] if j['dataset']==ds['id']),key=lambda j:j['id'])
        if not all((root/'chunks'/(j['id']+'.json')).exists() for j in jobs):print('INCOMPLETE '+ds['id']);fail=True;continue
        claim=root/'claims'/('merge_'+ds['id'])
        try:claim.mkdir()
        except FileExistsError:raise RuntimeError('Merge already claimed: '+ds['id'])
        try:
            output=root/'root'/('eppi_minus_'+ds['id']+'.root');report=root/'reports'/(ds['id']+'.json')
            counts=[read(root/'chunks'/(j['id']+'.json')) for j in jobs];entries=sum(c['entries'] for c in counts)
            if report.exists():
                if digest(output)!=read(report)['sha256']:raise RuntimeError('Published ROOT changed')
                print('Already merged '+ds['id']);continue
            paths=[root/'chunks'/(j['id']+'.root') for j in jobs]
            for p,c in zip(paths,counts):
                if digest(p)!=c['root_sha256']:raise RuntimeError('Chunk changed: '+str(p))
                run([build/'converter','--validate',p,str(c['entries'])])
            listing=root/'work'/('merge_'+ds['id']+'.list');listing.write_text(''.join(str(p)+'\n' for p in paths));tmp=output.with_suffix('.partial.root')
            run(['hadd','-f',tmp,'@'+str(listing)])
            run([build/'converter','--validate',tmp,str(entries)])
            os.replace(tmp,output);atomic_json(report,dict(file=str(output),entries=entries,batches=len(jobs),sha256=digest(output),zero_entries=entries==0))
            print(f'MERGED {output}: {entries} entries'+(' WARNING: EMPTY TREE' if entries==0 else ''),flush=True)
        finally:claim.rmdir()
    if fail:raise RuntimeError('Some datasets are incomplete; see status and batch logs')

def recover_conversion(a):
    """Reconvert completed text with a separately recorded converter repair."""
    root=a.campaign.resolve()
    if any((root/'hosts').iterdir()):raise RuntimeError('Wait for both worker commands to finish before conversion recovery')
    lock=root/'claims'/'conversion_recovery';lock.mkdir()
    try:
        build=verified(root);manifest=read(root/'manifest.json')
        recovery=root/'recoveries'/(time.strftime('%Y%m%d_%H%M%S')+'_'+str(os.getpid()))
        recovery.mkdir(parents=True)
        shutil.copytree(HERE/'src',recovery/'src')
        flags=shlex.split(subprocess.check_output(['root-config','--cflags','--libs'],text=True))
        converter=recovery/'converter'
        run(['g++','-O3','-std=c++17',recovery/'src/convert_exclusive_piminus_txt_to_root.cpp','-o',converter,*flags])
        env=os.environ.copy();env['PION_CONVERTER']=str(converter)
        run([sys.executable,'-m','unittest','discover','-s',HERE/'tests','-p','test_converter.py','-v'],env=env)
        atomic_json(recovery/'provenance.json',dict(converter_sha256=digest(converter),sources={str(p.relative_to(recovery)):digest(p) for p in (recovery/'src').rglob('*') if p.is_file()},original_software_unchanged=True))
        failures=[]
        for job in manifest['jobs']:
            jid=job['id'];done=root/'chunks'/(jid+'.json')
            if done.exists():continue
            claim=root/'claims'/jid;text=root/'work'/jid/'events.txt';log=root/'logs'/(jid+'.log')
            if not (claim/'failed.json').exists():
                failures.append(jid);print('NOT RECOVERED (no failed completion): '+jid,flush=True);continue
            summaries=re.findall(r'SUMMARY[^\n]*?rows=(\d+)',log.read_text(errors='replace')) if log.exists() else []
            if not summaries or not text.exists():
                failures.append(jid);print('NOT RECOVERED (processor did not finish): '+jid,flush=True);continue
            expected=int(summaries[-1]);tmp=root/'work'/jid/'recovered.root';started=time.time()
            result=subprocess.run([str(converter),str(text),str(tmp),str(build/'run_info.csv'),job['target'],str(job['run'])],capture_output=True,text=True)
            with log.open('a') as out:out.write('\nCONVERSION RECOVERY '+str(recovery)+'\n'+result.stdout+result.stderr)
            try:
                if result.returncode:raise RuntimeError(result.stderr.strip() or f'Converter exit {result.returncode}')
                entries=int(result.stdout.strip().splitlines()[-1])
                if entries!=expected:raise RuntimeError(f'Processor summary {expected} rows != converter {entries}')
                run([converter,'--validate',tmp,str(entries)])
                final=root/'chunks'/(jid+'.root');os.replace(tmp,final)
                atomic_json(done,dict(entries=entries,seconds=time.time()-started,host=socket.gethostname(),bytes=job['bytes'],root_sha256=digest(final),conversion_recovery=str(recovery)))
                (claim/'failed.json').unlink();(claim/'owner.json').unlink(missing_ok=True);claim.rmdir()
                print(f'RECOVERED {jid}: {entries} entries; conversion {time.time()-started:.1f}s',flush=True)
            except Exception as error:
                failures.append(jid);print(f'FAILED {jid}: {error}',flush=True)
        if failures:raise RuntimeError(f'{len(failures)} batches still need attention; saved text files retained')
        merge(a)
    finally:lock.rmdir()

def retry(a):
    root=a.campaign.resolve()
    if any((root/'hosts').iterdir()):raise RuntimeError('Host claims exist; confirm workers have stopped before recovering stale claims')
    for p in (root/'claims').iterdir():
        if (p/'failed.json').exists():
            shutil.rmtree(p);print('Released failed batch '+p.name)
    print('Failed batches released; run the worker command on each host again. Active/stale claims were not removed.')

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    for name in ['prepare','worker','status','merge','retry-failed','recover-conversion','monitor']:
        s=sub.add_parser(name);s.add_argument('--campaign',type=Path,required=True)
        if name=='prepare':
            s.add_argument('--repo',type=Path,required=True);s.add_argument('--mode',choices=['test','full'],default='test');s.add_argument('--batch-gib',type=float,default=4);s.add_argument('--run-map',type=Path);s.add_argument('--test-runs',type=Path)
        if name=='monitor':s.add_argument('--interval',type=int,default=30);s.add_argument('--once',action='store_true')
        if name=='worker':s.add_argument('--workers',type=int,choices=range(1,11),default=10);s.add_argument('--allow-other-host',action='store_true',help=argparse.SUPPRESS)
    a=p.parse_args()
    if a.command=='monitor' and a.interval<1:p.error('--interval must be positive')
    if a.command=='prepare' and a.batch_gib<=0:p.error('--batch-gib must be positive')
    {'prepare':prepare,'worker':worker,'status':lambda a:status(a.campaign.resolve()),'merge':merge,'retry-failed':retry,'recover-conversion':recover_conversion,'monitor':monitor}[a.command](a)
if __name__=='__main__':
    try:main()
    except Exception as e:print('ERROR: '+str(e),file=sys.stderr);sys.exit(1)
