"""Shared, validated production-response input routing and unchanged DATA cache."""
import hashlib,json,os,uuid
from pathlib import Path
import numpy as np
_state=None

def routed_files(files):
    global _state
    name=os.environ.get('DVCS_PARTICLE_RESPONSE_MANIFEST')
    if not name:return files
    if _state is None or _state[0]!=name:
        m=json.loads(Path(name).read_text());_state=(name,m,set())
    _,m,checked=_state;result=[]
    for name in files:
        key=os.path.abspath(name)
        if key in m['routes']:
            r=m['routes'][key];target=r['target']
            if key not in checked:
                for path,expected in [(key,r['source_identity']),(target,r['target_identity'])]:
                    s=Path(path).stat()
                    if (s.st_size,s.st_mtime_ns)!=(expected['size'],expected['mtime_ns']):raise ValueError('Response input changed: '+path)
                checked.add(key)
            result.append(target)
        elif key in m['preserved_files']:result.append(name)
        else:raise ValueError('Input missing from frozen production response manifest: '+key)
    return result

def cached_data_arrays(key,loader):
    """Only callers that identified DATA may use this cache; MC never shares it."""
    directory=os.environ.get('DVCS_PARTICLE_RESPONSE_DATA_CACHE')
    if not directory:return loader()
    digest=hashlib.sha256(json.dumps(key,sort_keys=True,default=str).encode()).hexdigest();folder=Path(directory)/digest
    marker=folder/'complete.json'
    if marker.is_file():
        meta=json.loads(marker.read_text());return {top:{name:np.load(folder/file,mmap_mode='r',allow_pickle=False) for name,file in vv.items()} for top,vv in meta.items()}
    result=loader();folder.mkdir(parents=True,exist_ok=True);meta={}
    for i,(top,values) in enumerate(result.items()):
        meta[top]={}
        for j,(name,a) in enumerate(values.items()):
            file=f'{i}_{j}.npy';tmp=folder/(file+'.'+uuid.uuid4().hex+'.tmp')
            with tmp.open('wb') as f:np.save(f,a,allow_pickle=False)
            tmp.replace(folder/file);meta[top][name]=file
    tmp=folder/('complete.'+uuid.uuid4().hex+'.tmp');tmp.write_text(json.dumps(meta));tmp.replace(marker);return result
