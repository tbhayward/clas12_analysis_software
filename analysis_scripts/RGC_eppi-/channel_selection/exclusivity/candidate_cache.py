"""One ROOT scan followed by small, memory-mapped diagnostic records."""
import hashlib
import inspect
import json
import os
from pathlib import Path
import shutil
import time
import numpy as np
from kinematics import BRANCHES, XB, TP, MN, MD, ME, MP, MPI, track, mass2, bin_indices

DTYPE=np.dtype([('mass','<f8'),('momentum','<f8'),('angle','<f8'),('bin','u1'),('detector','u1')])
CACHE_VERSION=1

def save_json(path,value):
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    tmp.replace(path)

def select(a):
    finite=np.logical_and.reduce([np.isfinite(a[k]) for k in BRANCHES])
    quality=(a['fiducial_status']==3)&(a['detector2']==1)&np.isin(a['detector1'],[1,2])&(a['Q2']>=1)&(a['W']>=2)&(a['y']<=.80)
    bins=bin_indices(a['x'],-a['tprime'])
    good=finite&quality&(bins>=0)
    return good,bins

def derive(a):
    """Compute only the three needed observables, preserving double precision."""
    e,p,pi=track(a,'e',ME),track(a,'p1',MP),track(a,'p2',MPI)
    pred=-e-p
    pred[:,0]+=a['beam_energy']+MN
    pred[:,3]+=np.sqrt(a['beam_energy']**2-ME**2)
    missing=pred-pi
    if len(missing) and np.max(np.abs(mass2(missing)-a['Mx2']))>1e-6:
        raise ValueError('Stored Mx2 does not match neutron-reference reconstruction')
    momentum=np.linalg.norm(missing[:,1:],axis=1)
    missing[:,0]+=MD-MN
    mass=mass2(missing)
    pmag=np.linalg.norm(pred[:,1:],axis=1)
    with np.errstate(divide='ignore',invalid='ignore'):
        theta=np.degrees(a['p2_theta']-np.arccos(np.clip(pred[:,3]/pmag,-1,1)))
    theta=np.where((pmag>1e-12)&(a['p2_p']>1e-12),theta,np.nan)
    return mass,momentum,theta

def build_cache(source,tree_name,out,options):
    import uproot
    source=Path(source)
    stat=source.stat()
    signature=dict(version=CACHE_VERSION,source=str(source.resolve()),size=stat.st_size,mtime_ns=stat.st_mtime_ns,tree=tree_name,max_entries=options['max_entries'],
                   selection_sha256=hashlib.sha256((inspect.getsource(select)+inspect.getsource(derive)+str(XB.tolist())+str(TP.tolist())+str((MN,MD,ME,MP,MPI))+str(BRANCHES)+inspect.getsource(track)+inspect.getsource(mass2)+inspect.getsource(bin_indices)).encode()).hexdigest())
    cache=out/'derived_cache'
    manifest=cache/'manifest.json'
    if manifest.exists() and not options.get('rebuild_cache',False):
        old=json.loads(manifest.read_text())
        if old['signature']==signature and all((cache/f).is_file() for f in old['chunks']):
            print(f'{out.name}: reused {old["combinations"]:,} cached candidates, no ROOT scan',flush=True)
            return cache,old
    temporary=out/f'derived_cache.build-{os.getpid()}'
    if temporary.exists():shutil.rmtree(temporary)
    temporary.mkdir()
    files=[]
    read=selected=0
    started=time.monotonic()
    save_json(out/'progress.json',dict(stage='ROOT cache',state='opening input',updated=time.time()))
    with uproot.open(source) as f:
        tree=f[tree_name]
        absent=set(BRANCHES)-set(tree.keys())
        if absent:raise ValueError(f'Missing branches {sorted(absent)}')
        total=min(tree.num_entries,options['max_entries'] or tree.num_entries)
        source_entries=tree.num_entries
        save_json(out/'progress.json',dict(stage='ROOT cache',state='reading first chunk',percent=0,updated=time.time()))
        for a in tree.iterate(BRANCHES,library='np',step_size=options['step_size'],entry_stop=total):
            read+=len(a['x'])
            good,bins=select(a)
            a={k:v[good] for k,v in a.items()}
            mass,momentum,angle=derive(a)
            records=np.empty(len(mass),dtype=DTYPE)
            records['mass'],records['momentum'],records['angle']=mass,momentum,angle
            records['bin'],records['detector']=bins[good],a['detector1']
            selected+=len(records)
            if len(records):
                name=f'chunk_{len(files):05d}.npy'
                np.save(temporary/name,records,allow_pickle=False)
                files.append(name)
            save_json(out/'progress.json',dict(stage='ROOT cache',state='scanning',percent=100*read/max(total,1),read=read,total=total,updated=time.time()))
            print(f'{out.name}: ROOT scan {100*read/max(total,1):.1f}%, {selected:,} candidates cached',flush=True)
    meta=dict(signature=signature,chunks=files,combinations=selected,source_entries=source_entries,read=read,partial=read<source_entries,elapsed_seconds=time.monotonic()-started,record_bytes=DTYPE.itemsize,angle_unit='degree')
    save_json(temporary/'manifest.json',meta)
    if cache.exists():shutil.rmtree(cache)
    temporary.rename(cache)
    return cache,meta

def iterate(cache,meta):
    for name in meta['chunks']:
        yield np.load(cache/name,mmap_mode='r',allow_pickle=False)
