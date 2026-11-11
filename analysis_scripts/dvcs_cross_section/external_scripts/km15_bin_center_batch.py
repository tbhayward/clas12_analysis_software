#!/usr/bin/env python3
"""Converged KM15 centering on the same allowed domain as bin volumes.
Stored Fbin is a denominator: model average over allowed bin / model(mean).
Use --output to preserve the input; without it the requested CSV is updated atomically.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,csv,hashlib,json,math,sys,tempfile,warnings
sys.dont_write_bytecode=True
os.environ["PYTHONDONTWRITEBYTECODE"]="1"
from multiprocessing import get_context
from pathlib import Path
import numpy as np
from kinematic_domain import Domain,box_from_row,energy_for_label,find_config
warnings.filterwarnings('ignore',category=RuntimeWarning)
_G=_TH=None
VERSION='allowed-domain-v3-pass1-divide'

def install_cff_cache(model):
    """Reuse fixed-KM15 CFF scalars at identical invariants, independently of phi.

    Initial reuses at a different phi are checked against the original method.
    The cached Fourier coefficients are phi independent by construction; the
    actual trigonometric and propagator factors are still evaluated per point.
    This optimizes repeated dispersion integrals, not model/grid interpolation.
    The model parameters remain fixed for the lifetime of a worker.
    """
    from collections import OrderedDict
    for name in ('ReH','ImH','ReE','ImE','ReHt','ImHt','ReEt','ImEt',
                 'cBH0unp','cBH1unp','cBH2unp','cINT0unp','cINT1unp','cINT2unp','cINT3unp',
                 'sINT1unp','sINT2unp','sINT3unp','cDVCS0unp','cDVCS1unp','sDVCS1unp'):
        original=getattr(model,name,None)
        if original is None:continue
        cache=OrderedDict();checks=[0]
        def cached(pt,*args,_original=original,_cache=cache,_checks=checks,**kwargs):
            try:
                key=(float(pt.xB),float(pt.Q2),float(pt.t),float(pt.xi),
                     pt.get('in2particle','p'),float(pt.in1energy),pt.get('in1charge',-1),pt.get('in1polarization',0),args,tuple(sorted(kwargs.items())))
                phi=float(pt.phi);hash(key)
            except (AttributeError,KeyError,TypeError):return _original(pt,*args,**kwargs)
            if key in _cache:
                value,first_phi,checked=_cache[key]
                if not checked and phi!=first_phi and _checks[0]<4:
                    reference=_original(pt,*args,**kwargs)
                    if float(reference)!=float(value):
                        raise ArithmeticError('KM15 CFF changed with phi; refusing cached evaluation')
                    _cache[key]=(value,first_phi,True);_checks[0]+=1
                _cache.move_to_end(key)
                return value
            value=_original(pt,*args,**kwargs);_cache[key]=(value,phi,False)
            if len(_cache)>4096:_cache.popitem(last=False)
            return value
        setattr(model,name,cached)
    return model

def _init_worker():
    global _G,_TH
    import gepard as g
    from gepard.fits import th_KM15
    import copy
    _G,_TH=g,copy.copy(th_KM15)
    # Boundary-aware nodes introduce many distinct Q² scales. Bound the
    # model's pure Wilson-kernel cache rather than growing it across all rows.
    from collections import OrderedDict
    class KernelCache(OrderedDict):
        def __getitem__(self,key):
            value=super().__getitem__(key);self.move_to_end(key);return value
        def __setitem__(self,key,value):
            super().__setitem__(key,value);self.move_to_end(key)
            if len(self)>256:self.popitem(last=False)
    kernels=KernelCache()
    for key,value in getattr(_TH,'wce',{}).items():kernels[key]=value
    _TH.wce=kernels
    if not os.environ.get('KM15_DISABLE_CFF_CACHE'):install_cff_cache(_TH)

def _km15(point):
    x,q,t,phi,e=point
    pt=_G.DataPoint(xB=x,t=-abs(t),Q2=q,phi=math.pi-math.radians(phi),observable='XS',
        frame='trento',process='ep2epgamma',exptype='fixed target',in1energy=e,
        in1charge=-1,in1polarization=0)
    pt.prepare();value=float(_TH.predict(pt))
    if not math.isfinite(value) or value<=0:raise ArithmeticError('Invalid KM15 value inside allowed domain')
    return value

def triple(s):
    try:
        p=[float(v.strip()) for v in str(s).strip().strip('"').strip('()').split(',')]
        return p[0],p[1] if len(p)>1 else 0.,p[2] if len(p)>2 else 0.
    except (ValueError,TypeError):return math.nan,math.nan,math.nan

def fmt(v,s=0.,y=0.):return f'({v:.12g}, {s:.12g}, {y:.12g})'

def integrate_task(job):
    i,label,center,box,phi,e,order,max_order,rtol,domain,volume=job
    last={}
    try:
        if not domain.contains(center):raise ArithmeticError('Measured mean is outside generated phase space')
        central=_km15(center);previous=None;cache={};last={}
        orders=sorted(set([order,order+1,*range(order+3,max_order+1,2),max_order]))
        for n in orders:
            points,weights=domain.quadrature(box,phi,e,n)
            if not len(points):raise ArithmeticError('No allowed quadrature points')
            vals=[]
            for point in points:
                if point not in cache:cache[point]=_km15(point)
                vals.append(cache[point])
            approximate_volume=float(weights.sum())
            average=float(np.dot(weights,np.asarray(vals))/approximate_volume)
            factor=average/central
            change=abs(factor/previous-1) if previous is not None else math.inf
            volume_error=abs(approximate_volume/volume-1)
            last=dict(row=i,label=label,order=n,factor=factor,model_at_mean=central,
                model_bin_average=average,relative_change=change,quadrature_volume_relative_error=volume_error,
                evaluations=len(cache),status='unconverged')
            if previous is not None and change<=rtol and volume_error<=rtol:
                last['status']='converged';return last
            previous=factor
        return last
    except Exception as exc:
        return dict(row=i,label=label,status='invalid',error=str(exc),**{k:v for k,v in last.items() if k not in ('row','label','status')})

def atomic_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=str(path.parent));os.close(fd)
    try:
        with open(name,'w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        os.replace(name,path)
    finally:
        if os.path.exists(name):os.unlink(name)

def recalculate(rows,domain,workers=10,order=3,max_order=11,rtol=.002,fbin=True):
    if not 1<=workers<=10:raise ValueError('workers must be 1–10')
    if not 2<=order<max_order<=15:raise ValueError('Require 2 <= order < max-order <= 15')
    if not 0<rtol<.1:raise ValueError('rtol must be positive and below 0.1')
    updated=[dict(r) for r in rows];jobs=[];diagnostics=[]
    for i,r in enumerate(rows):
        box=box_from_row(r);phi=(float(r['phimin']),float(r['phimax']))
        if phi[1]-phi[0]>=359:continue
        for label,means_label in [('10.6 GeV','10.6 GeV'),('10.2 GeV','Sp19 Inb')]:
            e=energy_for_label(float(label.split()[0]));volume,cubic=domain.volume(box,phi,e)
            updated[i]['bin_volume, '+label]=fmt(volume);updated[i]['cubic bin_volume, '+label]=fmt(cubic)
            diagnostics.append(dict(row=i,label=label,volume=volume,cubic_volume=cubic,allowed_fraction=volume/cubic,status='volume'))
            if not fbin:continue
            center=tuple(triple(r.get(p+', '+means_label,''))[0] for p in ('xBavg','Q2avg','t_abs_avg','phiavg'))+(e,)
            if not all(math.isfinite(v) for v in center):continue
            if volume<=0:raise ArithmeticError(f'Populated row {i} has zero allowed volume')
            jobs.append((i,label,center,box,phi,e,order,max_order,rtol,domain,volume))
    if jobs:
        print(f'[bincenter] {len(jobs)} populated energy-rows; {workers} workers; tolerance={rtol:g}',flush=True)
        if workers==1:
            _init_worker();results=map(integrate_task,jobs)
            for result in results:diagnostics.append(result)
        else:
            with get_context('spawn').Pool(workers,initializer=_init_worker) as pool:
                for count,result in enumerate(pool.imap_unordered(integrate_task,jobs,chunksize=1),1):
                    diagnostics.append(result)
                    if count%20==0 or count==len(jobs):print(f'[bincenter] {count}/{len(jobs)} energy-rows complete',flush=True)
        for d in diagnostics:
            if d['status']=='converged':updated[d['row']]['Fbin, '+d['label']]=fmt(d['factor'])
    return updated,diagnostics

def rescale_cross_sections(original,updated):
    """Exact frozen-calibration update: sigma_new/sigma_old = Fold/Fnew * Vold/Vnew."""
    for a,b in zip(original,updated):
        for col in a:
            if not ((col.startswith('cross sections,') or col.startswith('normed cross sections,')) and ', exp,' in col and col.endswith(('unpol','pos','neg'))):continue
            value,stat,syst=triple(a[col])
            if not math.isfinite(value):continue
            label='10.2 GeV' if 'Sp19 Inb' in col or '10.2 GeV' in col else '10.6 GeV'
            if b['Fbin, '+label]=='':
                b[col]=''
                continue
            oldf=triple(a['Fbin, '+label])[0];newf=triple(b['Fbin, '+label])[0]
            oldv=triple(a['bin_volume, '+label])[0];newv=triple(b['bin_volume, '+label])[0]
            if not all(math.isfinite(v) and v>0 for v in (oldf,newf,oldv,newv)):
                raise ArithmeticError(f'Cannot rescale populated cross section: {col}')
            ratio=oldf/newf*oldv/newv;b[col]=fmt(value*ratio,stat*abs(ratio),syst*abs(ratio))

CONVENTION_COLUMN='cross section correction convention'

def convert_input_convention(rows,requested=None):
    """Return a copy using Lee denominators; refuse to guess an unmarked CSV."""
    recorded={r.get(CONVENTION_COLUMN,'') for r in rows}
    if requested is not None and recorded != {''} and recorded != {requested}:
        raise ValueError('Requested correction convention conflicts with CSV markers')
    convention=requested or (next(iter(recorded)) if len(recorded)==1 else None)
    if convention not in ('legacy-multiply','pass1-divide'):
        raise ValueError('Unmarked/ambiguous CSV: specify --input-correction-convention legacy-multiply or pass1-divide')
    updated=[dict(r) for r in rows]
    for a,b in zip(rows,updated):
        if convention=='legacy-multiply':
            for col in a:
                if not ((col.startswith('cross sections,') or col.startswith('normed cross sections,')) and ', exp,' in col and col.endswith(('unpol','pos','neg'))):continue
                value,stat,syst=triple(a[col])
                if not math.isfinite(value):continue
                label='10.2 GeV' if 'Sp19 Inb' in col or '10.2 GeV' in col else '10.6 GeV'
                frad=triple(a['Frad, '+label])[0];fbin=triple(a['Fbin, '+label])[0]
                if not all(math.isfinite(v) and v>0 for v in (frad,fbin)):
                    raise ValueError(f'Invalid imported correction for {col}')
                ratio=1/(frad*fbin)**2
                b[col]=fmt(value*ratio,stat*ratio,syst*ratio)
        b[CONVENTION_COLUMN]='pass1-divide'
    return updated,convention

def write_diagnostics(path,diagnostics):
    fields=list(dict.fromkeys(k for d in diagnostics for k in d))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(diagnostics)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('csv_path',type=Path)
    ap.add_argument('--output',type=Path);ap.add_argument('--cuts-json',type=Path)
    ap.add_argument('--workers',type=int,default=10);ap.add_argument('--order',type=int,default=3)
    ap.add_argument('--max-order',type=int,default=11);ap.add_argument('--rtol',type=float,default=.002)
    ap.add_argument('--on-unconverged',choices=('error','omit'),default='error',help='error: preserve output on failure; omit: leave failed Fbin/XS energy-rows blank for a preview')
    ap.add_argument('--input-correction-convention',choices=('legacy-multiply','pass1-divide'))
    ap.add_argument('--volumes-only',action='store_true');ap.add_argument('--update-cross-sections',action='store_true')
    args=ap.parse_args();source=args.csv_path.resolve();target=(args.output or source).resolve()
    original_bytes=source.read_bytes()
    import io
    rows=list(csv.DictReader(io.StringIO(original_bytes.decode('utf-8-sig'))))
    if not rows:raise SystemExit('CSV has no data rows')
    if args.update_cross_sections:
        try:rows,input_convention=convert_input_convention(rows,args.input_correction_convention)
        except ValueError as e:ap.error(str(e))
    config=find_config(source,args.cuts_json);domain=Domain.from_json(config)
    print(f'[bincenter] Domain configuration: {config}',flush=True)
    updated,diagnostics=recalculate(rows,domain,args.workers,args.order,args.max_order,args.rtol,not args.volumes_only)
    target.parent.mkdir(parents=True,exist_ok=True);diagpath=Path(str(target)+'.numerical_diagnostics.csv');write_diagnostics(diagpath,diagnostics)
    failed=[d for d in diagnostics if d['status'] not in ('converged','volume')]
    if failed and args.on_unconverged=='error':raise SystemExit(f'{len(failed)} energy-rows failed convergence/model checks; CSV not replaced. Inspect {diagpath}')
    for d in failed:updated[d['row']]['Fbin, '+d['label']]=''
    if failed:print(f'[bincenter] Omitting {len(failed)} unresolved Fbin/XS energy-rows from this preview.',flush=True)
    if args.update_cross_sections:rescale_cross_sections(rows,updated)
    meta=dict(version=VERSION,input=str(source),input_sha256=hashlib.sha256(original_bytes).hexdigest(),
        config=str(config),order=args.order,max_order=args.max_order,rtol=args.rtol,
        reference='imported' if args.volumes_only else 'measured_means',cross_sections_updated=args.update_cross_sections,
        normalization='held fixed',systematics='not reevaluated',on_unconverged=args.on_unconverged,omitted_energy_rows=[dict(row=d['row'],bin_index=rows[d['row']].get('bin index',''),label=d['label'],status=d['status']) for d in failed])
    if source.read_bytes()!=original_bytes:raise SystemExit('Input CSV changed during calculation; output not replaced')
    atomic_csv(target,updated);Path(str(target)+'.numerical.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(f'[bincenter] Saved {target}; original preserved: {source!=target}',flush=True)
    if not args.update_cross_sections:
        print('[bincenter] Cross-section columns unchanged; use --update-cross-sections or run_numerical_corrections.py for updated results.',flush=True)

if __name__=='__main__':main()
