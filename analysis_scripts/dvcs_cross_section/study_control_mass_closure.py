#!/usr/bin/env python3
"""Check saved versus rebuilt pi0 mass on existing per-event control exports."""
from __future__ import annotations
import argparse,html,json
from pathlib import Path
import numpy as np
import study_particle_response as base
import study_unordered_control_response as control
from study_particle_common_cut_stability import export_digest


def rebuilt_mass(a):
    # Stable 1-cos(opening) expression, avoiding cancellation for collinear photons.
    p1,t1,f1=a[:,20],a[:,21],a[:,22];p2,t2,f2=a[:,24],a[:,25],a[:,26]
    valid=np.all(np.isfinite(a[:,[20,21,22,24,25,26,28]]),axis=1)&(p1>0)&(p2>0)
    term=np.sin((t1-t2)/2)**2+np.sin(t1)*np.sin(t2)*np.sin((f1-f2)/2)**2
    mass=np.sqrt(np.maximum(0.,4*p1*p2*term));return mass,valid


def metrics(a,mask):
    rebuilt,valid=rebuilt_mass(a);use=mask&valid;saved=a[use,28];rebuilt=rebuilt[use];delta=1000*(rebuilt-saved)
    values=dict(selected=int(mask.sum()),resolved=int(use.sum()),unresolved=int((mask&~valid).sum()))
    if not len(delta):return dict(values,median_delta_MeV=float('nan'),p01_delta_MeV=float('nan'),p99_delta_MeV=float('nan'),rms_delta_MeV=float('nan'),max_abs_delta_MeV=float('nan'),fraction_above_0p1_MeV=float('nan'),fraction_above_1_MeV=float('nan'),saved_core_width_MeV=float('nan'),rebuilt_core_width_MeV=float('nan')),saved,rebuilt,delta
    q=np.quantile(delta,[.01,.5,.99]);sq=np.quantile(saved,[.16,.84]);rq=np.quantile(rebuilt,[.16,.84])
    return dict(values,median_delta_MeV=q[1],p01_delta_MeV=q[0],p99_delta_MeV=q[2],rms_delta_MeV=np.sqrt(np.mean(delta**2)),max_abs_delta_MeV=np.max(np.abs(delta)),fraction_above_0p1_MeV=np.mean(np.abs(delta)>.1),fraction_above_1_MeV=np.mean(np.abs(delta)>1),saved_core_width_MeV=500*(sq[1]-sq[0]),rebuilt_core_width_MeV=500*(rq[1]-rq[0])),saved,rebuilt,delta


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--control-study',required=True,type=Path);p.add_argument('--output',required=True,type=Path);args=p.parse_args(argv);origin=args.control_study.resolve();out=args.output.resolve()
    if out.exists():p.error('Output exists; choose a fresh directory')
    em=json.loads((origin/'export_manifest.json').read_text());spec=json.loads((origin/'spec.json').read_text());records=em['records'];expected={(period,scenario) for period in spec['periods'] for scenario in ('DATA','nominal','frozen')}
    if len(records)!=len(expected) or {(r['period'],r['scenario']) for r in records}!=expected:p.error('Incomplete control export inventory')
    for r in records:
        shape=r['shape'];path=Path(r['path'])
        if len(shape)!=2 or shape[1]!=29 or path.stat().st_size!=int(np.prod(shape))*8 or export_digest(path)!=r['sha256']:p.error('Export missing or changed: '+str(path))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out.mkdir(parents=True);rows=[];curves={};links=[]
    for r in records:
        period,scenario=r['period'],r['scenario'];a=np.memmap(r['path'],dtype='<f8',mode='r',shape=tuple(r['shape'])) if r['shape'][0] else np.empty((0,29))
        for group in control.ACTIVE_GROUPS:
            for scope,windows in [('global_only',[]),('common_windows',control.blocks(spec,period,group))]:
                mask=control.selection(a,group,windows);m,saved,rebuilt,delta=metrics(a,mask);rows.append(dict(period=period,group=control.GROUPS[group],sample=scenario,scope=scope,**m))
                if scope=='common_windows':
                    edges=np.linspace(.10,.17,201);n=np.histogram(saved,edges)[0];b=np.histogram(rebuilt,edges)[0];curves[period,group,scenario]=(n,b,edges,delta)
        print('[closure] '+period+'/'+scenario+'; existing arrays only',flush=True)
    for period in spec['periods']:
        for group in control.ACTIVE_GROUPS:
            fig,axes=plt.subplots(2,3,figsize=(14,7),layout='constrained')
            for j,scenario in enumerate(('DATA','nominal','frozen')):
                n,b,edges,delta=curves[period,group,scenario];ax=axes[0,j]
                for label,h,color in [('Saved/reported mass',n,'#205cab'),('Rebuilt from exported photons',b,'#d45c18')]:
                    if h.sum()>0:ax.stairs(h/h.sum()/np.diff(edges),edges,label=label,color=color)
                ax.set(title=scenario,xlabel='Two-photon mass (GeV)',ylabel='Unit-area density')
                if n.sum()+b.sum()>0:ax.legend(fontsize=7)
                ax=axes[1,j];limit=min(3.,max(.05,float(np.quantile(np.abs(delta),.999))*1.1)) if len(delta) else 3.
                ax.hist(delta,bins=np.linspace(-limit,limit,121),histtype='step',color='#d45c18');ax.axvline(0,color='black',lw=.7);ax.set(xlabel='Rebuilt − saved mass (MeV)',ylabel='Events',xlim=(-limit,limit));ax.grid(alpha=.2)
                if len(delta):ax.text(.02,.96,f'N={len(delta):,}\nmedian={np.median(delta):+.4g} MeV\noutside view={np.mean(np.abs(delta)>limit):.2%}',transform=ax.transAxes,va='top',fontsize=8)
            fig.suptitle(period+' · '+control.GROUPS[group]+' · common windows, matched resolved events\nDATA/nominal: zero-response reconstruction check; frozen: reported-mass consistency only');folder=out/'plots'/period;folder.mkdir(parents=True,exist_ok=True);path=folder/(control.GROUPS[group]+'.png');fig.savefig(path,dpi=110);plt.close(fig);links.append(path.relative_to(out).as_posix())
    base.table(out/'mass_closure.csv',rows);base.dump(out/'manifest.json',dict(control_export_sha256=base.digest(origin/'export_manifest.json'),script_sha256=base.digest(Path(__file__)),source_study=str(origin)))
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Control mass reconstruction closure</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}</style><h1>Saved versus rebuilt pi0 mass</h1><p>No ROOT reads, fits or new response modifications. Uses the existing per-event photon coordinates and saved/reported mass. DATA and nominal MC test zero-response reconstruction on the real samples; frozen MC tests agreement between the reported modified mass and its exported modified photons. Every overlay uses the exact same resolved events on each side. Unresolved events are counted separately, not silently treated as closure.</p><p>The stored mass was an input to photon recovery, so this is numerical reconstruction consistency, not independent validation of a physical mass model. Closure must be small compared with the peak shift/width under investigation; no acceptance threshold is assumed automatically. Both global and simultaneous common-window summaries are in the CSV; plots show common windows. No FD-proton group containing FT photons is plotted. Existing source processing/mass-window truncation remains.</p><p><a href="mass_closure.csv">Mass differences, quantiles and widths</a></p>'+''.join('<p><a href="'+x+'">'+html.escape(x)+'</a></p>' for x in links));print('[done] '+str(out/'REPORT.html'),flush=True)

if __name__=='__main__':main()
