"""Full-statistics publication histograms, streamed once per MC file/variant."""
import concurrent.futures,hashlib,json,os,uuid
from pathlib import Path
import numpy as np
import uproot
import plot_exclusivity_data_dvcs_pi0_mc as production
import study_particle_response as base

def stamp(path):
    s=Path(path).stat();return dict(path=str(path),size=s.st_size,mtime_ns=s.st_mtime_ns)
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True).encode()).hexdigest()
def mix(x):
    with np.errstate(over='ignore'):
        x=x+np.uint64(0x9e3779b97f4a7c15);x=(x^(x>>np.uint64(30)))*np.uint64(0xbf58476d1ce4e5b9);x=(x^(x>>np.uint64(27)))*np.uint64(0x94d049bb133111eb);return x^(x>>np.uint64(31))
def halves(path,seed,start,n,direct):
    key=1469598103934665603
    for c in path.encode():key=((key^c)*1099511628211)&((1<<64)-1)
    k=mix(np.arange(start,start+n,dtype=np.uint64)^np.uint64(key^seed));k=(k&np.uint64(0xfffffffffffffffe))|np.uint64(not direct)
    return ((mix(k^np.uint64(0xb234fe3adbf120dd))>>np.uint64(40))&np.uint64(1)).astype(np.int8)
def window_mask(values,window,omit=None):
    keep=np.ones(len(next(iter(values.values()))),dtype=bool)
    for name,w in window.items():
        if name==omit:continue
        v=values[name];keep&=np.isfinite(v)&(v<=float(w['cut_high']))
        if w.get('mode')!='upper_quantile':keep&=v>=float(w['cut_low'])
    return keep

def collect(manifest,blocks,variables,mass,out,workers=10):
    out=Path(out);jobs=[]
    for source,r in manifest['routes'].items():
        roles=[v['sample'] for v in r['roles'] if v['sample'] in ('mc','background','pi0_mc')]
        if not roles:continue
        sample=roles[0]
        for scenario,path in enumerate((source,r['target'])):jobs.append(dict(source=source,path=path,period=r['period'],sample=sample,scenario=scenario))
    recipe=dict(inputs=[dict(j,identity=stamp(j['path'])) for j in jobs],windows=blocks,globals=manifest['global_cuts'],variables=variables,mass=mass,seed=manifest['model']['seed'],collector=base.digest(Path(__file__)),global_reader=base.digest(Path(production.__file__)))
    key=digest(recipe);meta=out/'full_statistics_manifest.json';cache=out/'full_statistics_histograms.npz';periods=manifest['model']['periods'];samples=('mc','background','pi0_mc')
    if meta.is_file() and cache.is_file() and json.loads(meta.read_text()).get('signature')==key:
        with np.load(cache,allow_pickle=False) as z:result={p:({s:z[p+'_'+s].copy() for s in samples},{s:z[p+'_selected_'+s].copy() for s in samples}) for p in periods}
        print('[response note] Full-statistics MC histograms verified and reused; no ROOT reread',flush=True);return result
    result={}
    for period in periods:
        h={};selected={}
        for sample in samples:
            vv=variables+([mass] if sample=='pi0_mc' else []);h[sample]=np.zeros((2,2,3,len(vv),max(v['bins']+2 for v in vv)),dtype=np.int64);selected[sample]=np.zeros((2,2,3),dtype=np.int64)
        result[period]=(h,selected)
    def one(job):
        control=job['sample']=='pi0_mc';channel='epgammagamma' if control else 'epgamma';vv=variables+([mass] if control else []);hist=np.zeros((2,3,len(vv),max(v['bins']+2 for v in vv)),dtype=np.int64);counts=np.zeros((2,3),dtype=np.int64);executor=uproot.source.futures.TrivialExecutor();offset=0
        physical=['x','Q2','t1','detector1','detector2','e_p','e_theta','e_phi','p1_p','p1_theta','p1_phi','p2_p','p2_theta','p2_phi','open_angle_ep2','pTmiss','runnum'];observables=['theta_pi0_pi0' if control and n=='theta_gamma_gamma' else n for n in base.OBS];extra=['detector_gamma1','Mh_gammagamma'] if control else []
        with uproot.open(job['path'],decompression_executor=executor,interpretation_executor=executor,num_workers=1) as f:
            tree=f['PhysicsEvents'];direct='Delta_phi' in tree;fields=[n for n in dict.fromkeys(physical+observables+extra) if n in tree]
            for chunk in tree.iterate(expressions=fields,step_size=65536,library='np'):
                a=dict(chunk) if isinstance(chunk,dict) else {n:chunk[n] for n in chunk.dtype.names};n=len(next(iter(a.values())));a.setdefault('runnum',np.zeros(n))
                if control:a['detector2']=a['detector_gamma1']
                if not direct:a['Delta_phi']=np.abs((a['p1_phi']-a['p2_phi']+np.pi)%(2*np.pi)-np.pi)
                values={name:a['theta_pi0_pi0' if control and name=='theta_gamma_gamma' else name] for name in base.OBS}
                if control:values['m_gamma_gamma']=a['Mh_gammagamma']
                valid=production.global_cut_mask_uproot(a,job['period'],a['detector1'],a['detector2']);half=halves(job['source'],manifest['model']['seed'],offset,n,direct);offset+=n
                for g,(det1,det2) in enumerate(((1,1),(2,1),(2,0))):
                    pop=valid&(a['detector1']==det1)&(a['detector2']==det2);block=blocks[job['period']][channel][g]['data'];full=pop&window_mask(values,block)
                    for hh in (0,1):
                        which=half==hh;counts[hh,g]+=np.count_nonzero(full&which)
                        for j,v in enumerate(vv):
                            keep=pop&which&window_mask(values,block,v['name']);x=values[v['name']][keep];x=x[np.isfinite(x)];edges=base.diag.variable_edges(v);hist[hh,g,j,:len(edges)+1]+=np.bincount(np.searchsorted(edges,x,side='right'),minlength=len(edges)+1)
            if offset!=tree.num_entries:raise ValueError('Incomplete full-statistics publication scan')
        return job,hist,counts,offset
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(10,workers)) as pool:
        futures=[pool.submit(one,j) for j in jobs]
        for k,f in enumerate(concurrent.futures.as_completed(futures),1):
            job,h,c,n=f.result();hh,ss=result[job['period']];hh[job['sample']][job['scenario']]+=h;ss[job['sample']][job['scenario']]+=c;print(f"[response note full MC] {k}/{len(jobs)} files; {job['period']}/{job['sample']}/{'nominal' if job['scenario']==0 else 'modified'}: {n:,} events",flush=True)
    payload={}
    for p,(h,s) in result.items():
        for sample in samples:payload[p+'_'+sample]=h[sample];payload[p+'_selected_'+sample]=s[sample]
    tmp=cache.with_name(cache.name+'.'+uuid.uuid4().hex+'.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,**payload)
    tmp.replace(cache);meta.write_text(json.dumps(dict(signature=key,recipe=recipe,statistics='all available reconstructed MC events; no subsampling'),indent=2)+'\n');return result
