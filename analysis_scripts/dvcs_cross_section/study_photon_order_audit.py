#!/usr/bin/env python3
"""Audit photon role ordering and control topology using existing frozen exports."""
from __future__ import annotations
import argparse,csv,html,json
from pathlib import Path
import numpy as np
import study_particle_response as base
import study_response_cut_redetermination as cuts
import study_particle_selection_closure as closure
from study_particle_common_cut_stability import export_digest

TOPS=['FD_FD','CD_FD','CD_FT','FD_FT (not in current grouping)']

def summarize(a,period,sample,scope,keep):
    tags=np.isin(a[:,23],[0,1])&np.isin(a[:,27],[0,1]);resolved=tags&np.isfinite(a[:,20])&np.isfinite(a[:,24])&(a[:,20]>0)&(a[:,24]>0)
    valid=keep&resolved;hard1=a[:,20]>a[:,24];ties=a[:,20]==a[:,24]
    rows=[];migration=[]
    for d1,d2,label in [(1,1,'FD_FD photons'),(1,0,'FD_FT photons'),(0,1,'FT_FD photons'),(0,0,'FT_FT photons'),(-1,-1,'all resolved')]:
        m=valid if d1<0 else valid&(a[:,23]==d1)&(a[:,27]==d2);n=int(m.sum());n1=int((m&hard1).sum());nt=int((m&ties).sum())
        rows.append(dict(period=period,sample=sample,scope=scope,photon_pair=label,events=n,photon1_harder=n1,photon2_harder=n-n1-nt,equal_energy=nt,photon1_harder_fraction=n1/n if n else float('nan'),scope_events=int(keep.sum()),unresolved_events=int((keep&~resolved).sum())))
    # Photon coordinates and detector tags are kept together. Do not relabel source rows.
    hard_tag=np.where(hard1,a[:,23],a[:,27]);fd_proton=a[:,0]==0
    leading_top=np.where(fd_proton,np.where(hard_tag==1,0,3),np.where(hard_tag==1,1,2))
    for old in range(3):
        m=valid&~ties&(a[:,0]==old);n=int(m.sum())
        for new in range(4):
            count=int((m&(leading_top==new)).sum());migration.append(dict(period=period,sample=sample,scope=scope,saved_topology=TOPS[old],harder_photon_topology=TOPS[new],events=count,row_events=n,row_fraction=count/n if n else float('nan')))
    return rows,migration


def plot(period,rows,migrations,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    samples=['DATA','nominal_MC','frozen_MC'];fig,axes=plt.subplots(2,3,figsize=(15,8),layout='constrained')
    for j,sample in enumerate(samples):
        rr=[r for r in rows if r['period']==period and r['sample']==sample and r['scope']=='global'];ax=axes[0,j];labels=['FD_FD photons','FD_FT photons','FT_FD photons','FT_FT photons'];values=[next(r['photon1_harder_fraction'] for r in rr if r['photon_pair']==label) for label in labels];counts=[next(r['events'] for r in rr if r['photon_pair']==label) for label in labels]
        ax.bar(range(4),np.asarray(values)*100,color='#205cab' if sample=='DATA' else '#d45c18');ax.set(xticks=range(4),xticklabels=['FD, FD','FD, FT','FT, FD','FT, FT'],ylim=(0,110),ylabel='Photon 1 harder (%)',title=sample+' · photon detector pair');ax.grid(axis='y',alpha=.2)
        for k,(value,count) in enumerate(zip(values,counts)):ax.text(k,103,'N='+str(count),ha='center',fontsize=8)
        rr=[r for r in migrations if r['period']==period and r['sample']==sample and r['scope']=='global'];matrix=np.array([[next(r['row_fraction'] for r in rr if r['saved_topology']==old and r['harder_photon_topology']==new) for new in TOPS] for old in TOPS[:3]])*100
        ax=axes[1,j];im=ax.imshow(np.ma.masked_invalid(matrix),vmin=0,vmax=100,cmap='Blues',aspect='auto');ax.set(xticks=range(4),xticklabels=['FD–FD','CD–FD','CD–FT','FD–FT'],yticks=range(3),yticklabels=TOPS[:3],xlabel='Topology if harder photon defines label',ylabel='Saved topology',title='Descriptive regrouping: row percentages')
        for index,value in np.ndenumerate(matrix):
            if np.isfinite(value):ax.text(index[1],index[0],f'{value:.1f}%',ha='center',va='center',color='white' if value>65 else 'black',fontsize=8)
    fig.suptitle(period+' · control photon ordering audit\nGlobally retained events only; no cut reevaluation or rate normalization');folder=out/'plots';folder.mkdir(exist_ok=True);path=folder/(period+'.png');fig.savefig(path,dpi=120);plt.close(fig);return path.relative_to(out).as_posix()


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--frozen-study',required=True,type=Path);p.add_argument('--output',required=True,type=Path);args=p.parse_args(argv);study=args.frozen_study.resolve();out=args.output.resolve()
    if out.exists():p.error('Output exists; choose a new directory')
    em=json.loads((study/'export_manifest.json').read_text());spec=json.loads((study/'spec.json').read_text());records=[r for r in em['records'] if r['sample'] in ('pi0_data','pi0_mc')];expected={(period,scenario,sample) for period in spec['periods'] for scenario,sample in [('shared','pi0_data'),('nominal','pi0_mc'),('candidate','pi0_mc')]}
    if {(r['period'],r['scenario'],r['sample']) for r in records}!=expected or len(records)!=len(expected):p.error('Requires all frozen DATA, nominal MC and candidate MC control exports')
    for r in records:
        path=Path(r['path']);shape=r['shape']
        if len(shape)!=2 or shape[1]!=28 or path.stat().st_size!=int(np.prod(shape))*8 or export_digest(path)!=r['sha256']:p.error('Invalid or changed export: '+str(path))
    out.mkdir(parents=True);rows=[];migrations=[]
    for r in records:
        period=r['period'];sample={'shared':'DATA','nominal':'nominal_MC','candidate':'frozen_MC'}[r['scenario']];a=np.memmap(r['path'],dtype='<f8',mode='r',shape=tuple(r['shape']));cfg=next(c for c in cuts.production.PERIODS if c.key==period);global_keep=a[:,8]>.5;selected=np.zeros(len(a),dtype=bool)
        for top in range(3):
            block=spec['exclusivity'][cuts.production.cut_json_key('eppi0',cfg,cuts.production.TOPOLOGIES[top])];selected|=closure.mask(a,top,block,'data')
        for scope,keep in [('global',global_keep),('common_DATA_windows',selected)]:
            rr,mm=summarize(a,period,sample,scope,keep);rows+=rr;migrations+=mm
        print('[audit] '+period+'/'+sample+'; no ROOT reads',flush=True)
    base.table(out/'ordering.csv',rows);base.table(out/'topology_regrouping.csv',migrations);links=[plot(period,rows,migrations,out) for period in spec['periods']]
    base.dump(out/'manifest.json',dict(frozen_export_manifest_sha256=base.digest(study/'export_manifest.json'),fingerprint=em['fingerprint'],source_sha256=base.digest(Path(__file__)),read_only=True))
    (out/'REPORT.html').write_text('<!doctype html><meta charset="utf-8"><title>Photon order audit</title><style>body{font:17px system-ui;max-width:1100px;margin:35px auto}</style><h1>Photon order and control topology</h1><p>DATA, unmodified MC and frozen-response MC use the checked 28-column exports. No fits, ROOT reads or production changes. Detector pairs refer to each photon\'s own tag; the tags travel with their coordinates. Unknown reconstructed photon pairs remain counted as unresolved and do not enter ordering/regrouping.</p><p>Top row: which photon is harder, conditional on detector pair. Bottom row: how saved topology groups would change if the harder photon defined the control label in both DATA and MC. All matrices are row-normalized, not cross-section or absolute-yield comparisons. FD–FT is shown explicitly although not part of the current three-group scheme.</p><p>These are descriptive regroupings of already retained events. Global cuts and exclusivity windows were not reevaluated under a new grouping; events lost by the old selection cannot be recovered here. This is not a corrected selection or new contamination calculation. Common-DATA-window results in the CSV apply each event\'s saved topology windows before regrouping. Resolving mixed-detector grouping requires choosing an explicitly consistent control definition, not simply swapping tags alone.</p><p><a href="ordering.csv">Counts and unresolved coverage</a> · <a href="topology_regrouping.csv">Regrouping counts and fractions</a></p>'+''.join('<p><a href="'+x+'">'+html.escape(x)+'</a></p>' for x in links))
    print('[done] '+str(out/'REPORT.html'),flush=True)

if __name__=='__main__':main()
