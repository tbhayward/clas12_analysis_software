#!/usr/bin/env python3
"""Paired DVCS MC proton response: event-level and MC::RecMatch particle diagnostics.
No MC::Event.weight; candidate-level, not final DVCS selection. All angles in degrees.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uproot

M=0.9382720813

def statsplot(ax,x,y,edges,label):
    centers=.5*(edges[:-1]+edges[1:]); ix=np.searchsorted(edges,x,side='right')-1
    ns=[]; means=[]; med=[]; std=[]
    for j in range(len(centers)):
        v=y[(ix==j)&np.isfinite(y)]; ns.append(len(v)); means.append(np.mean(v) if len(v)>20 else np.nan); med.append(np.median(v) if len(v)>20 else np.nan); std.append(np.std(v) if len(v)>20 else np.nan)
    #endfor
    ax.plot(centers,means,'o-',ms=3,label='Mean');ax.plot(centers,med,'s-',ms=3,label='Median');ax.axhline(0,color='black',lw=1);ax.set(xlabel=label,ylabel='Reconstructed minus generated');ax.grid(alpha=.2);ax.legend()
    return np.asarray(ns),np.asarray(means),np.asarray(std)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',default='/volatile/clas12/thayward/dvcs_reprocessing_20261022/production/fa18_inb/dvcs_rec/50nA')
    p.add_argument('--output',default='/u/home/thayward/clas12_analysis_software/analysis_scripts/misc/diagnostic_dvcs')
    p.add_argument('--max-events',type=int,default=10000000)
    p.add_argument('--chunk',type=int,default=200000)
    p.add_argument('--max-particle-entries',type=int,default=50000000,help='Hard safety limit for particle rows scanned')
    a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    root=Path(a.input);folder=root if root.is_dir() else root.parent
    ef=root if root.is_file() else folder/'all.mc_events.root';pf=folder/'all.mc_particles.root'
    if not ef.exists() or not pf.exists():raise FileNotFoundError(f'Need both {ef} and {pf}')
    print('Event ROOT:',ef,'\nParticle ROOT:',pf,flush=True)
    with uproot.open(ef) as f:
        t=f['MCEvents'];fields=['event_sequence','gen_kinematics_valid','rec_kinematics_valid','rec_final_state','gen_t1','rec_t1','gen_p1_p','rec_p1_p','gen_p1_theta','rec_p1_theta','rec_detector1']
        missing=[z for z in fields if z not in t.keys()]
        if missing:raise ValueError(f'Missing event branches {missing}')
        d={k:[] for k in ['g','r','gp','rp','gt','rt','det']};maxn=min(a.max_events,t.num_entries)
        # Use the actual last selected event identifier, not an assumed 1-based entry index.
        selected_last_sequence=int(t['event_sequence'].array(entry_start=maxn-1,entry_stop=maxn,library='np')[0])
        selected_first_sequence=int(t['event_sequence'].array(entry_start=0,entry_stop=1,library='np')[0])
        print(f'Selected event_sequence range: {selected_first_sequence:,} to {selected_last_sequence:,}',flush=True)
        for start in range(0,maxn,a.chunk):
            b=t.arrays(fields,entry_start=start,entry_stop=min(start+a.chunk,maxn),library='np');v=(b['gen_kinematics_valid']>0)&(b['rec_kinematics_valid']>0)&(b['rec_final_state']>0)&(b['gen_t1']<0)&(b['rec_t1']<0)
            for k,field in [('g','gen_t1'),('r','rec_t1'),('gp','gen_p1_p'),('rp','rec_p1_p'),('gt','gen_p1_theta'),('rt','rec_p1_theta'),('det','rec_detector1')]:d[k].append(b[field][v])
            #endfor
            print(f'Events {min(start+a.chunk,maxn):,}/{maxn:,}',flush=True)
        #endfor
    #endwith
    d={k:np.concatenate(v) for k,v in d.items()};g=-d['g'];r=-d['r'];delta=r-g
    fig,axs=plt.subplots(2,2,figsize=(13,9));statsplot(axs[0,0],g,delta,np.linspace(.10,.60,51),'Generated |t| (GeV²)');axs[0,0].set_ylabel('Δ|t| (GeV²)')
    statsplot(axs[0,1],d['gp'],d['rp']-d['gp'],np.linspace(0,1.5,51),'Generated proton p (GeV)');axs[0,1].set_ylabel('Δp (GeV)')
    statsplot(axs[1,0],np.degrees(d['gt']),np.degrees(d['rt']-d['gt']),np.linspace(20,160,57),'Generated proton θ (deg)');axs[1,0].set_ylabel('Δθ (deg)')
    axs[1,1].hist(delta[(g>=.11)&(g<.25)],bins=120,range=(-.15,.15),histtype='step');axs[1,1].set(xlabel='Δ|t| (GeV²)',ylabel='Events',title='Generated 0.11 < |t| < 0.25 (GeV²)');axs[1,1].axvline(0,color='black')
    fig.suptitle('DVCS proton response: event-level reconstructed candidates');fig.tight_layout();fig.savefig(out/'v8_01_proton_response.png',dpi=170);plt.close(fig)
    # Detector1 convention must be verified against processing classes; display values without assuming codes.
    fig,ax=plt.subplots(figsize=(9,5))
    for code in np.unique(d['det']):
        v=(d['det']==code)&(g>=.11)&(g<.6)
        if np.sum(v)>100:statsplot(ax,g[v],delta[v],np.linspace(.11,.6,35),'Generated |t| (GeV²)');break
        #endif
    #endfor
    ax.clear()
    for code in np.unique(d['det']):
        v=(d['det']==code)&(g>=.11)&(g<.6)
        if np.sum(v)<100:continue
        #endif
        bins=np.linspace(.11,.6,35);ix=np.searchsorted(bins,g[v],side='right')-1;cent=.5*(bins[:-1]+bins[1:]);mean=[np.mean(delta[v][ix==j]) if np.sum(ix==j)>30 else np.nan for j in range(len(cent))];ax.plot(cent,mean,'o-',ms=3,label=f'detector1={code} (N={np.sum(v):,})')
    #endfor
    ax.axhline(0,color='black');ax.set(xlabel='Generated |t| (GeV²)',ylabel='Mean Δ|t| (GeV²)',title='Response by reconstructed proton detector1 code');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'v8_02_detector_response.png',dpi=170);plt.close(fig)
    # Particle tree is one entry per MC particle, not one per event. Restrict by event_sequence to the event cap.
    with uproot.open(pf) as f:
        trees={k.split(';')[0]:f[k] for k in f.keys() if isinstance(f[k],uproot.behaviors.TTree.TTree)}
        print('Particle ROOT trees:',{k:len(v.keys()) for k,v in trees.items()})
        if not trees:raise RuntimeError('No TTree in particle ROOT')
        t=next((v for v in trees.values() if 'truth_pid' in v.keys()),None)
        if t is None:raise RuntimeError(f'No truth_pid branch; available {[(k,list(v.keys())) for k,v in trees.items()]}')
        keys=['event_sequence','truth_pid','rec_pindex','rec_pid','gen_p','rec_p','gen_theta','rec_theta','gen_t_proton','rec_t_proton_raw','rec_t_proton_corrected','fitted_selected','rec_detector']
        missing=[k for k in keys if k not in t.keys()]
        if missing:raise RuntimeError(f'Missing particle branches: {missing}')
        records={k:[] for k in keys}; n=0; stopped_at_boundary=False; previous_last=None
        if selected_last_sequence < selected_first_sequence:
            raise RuntimeError('Event identifiers are not ordered at the endpoints; cannot safely select particle entries.')
        # Particle rows are ordered by event_sequence in the intended production format.
        # Check ordering as we go and stop at the first row beyond the selected event range.
        for b in t.iterate(keys,step_size=a.chunk,entry_stop=min(t.num_entries,a.max_particle_entries),library='np'):
            seq=b['event_sequence']
            if len(seq)==0:continue
            #endif
            if (previous_last is not None and seq[0]<previous_last) or np.any(seq[1:]<seq[:-1]):
                raise RuntimeError('Particle event_sequence is not monotonic. Cannot use early termination safely; inspect ROOT event mapping.')
            #endif
            previous_last=seq[-1]
            eligible=seq<=selected_last_sequence
            sel=(b['truth_pid']==2212)&(seq>=selected_first_sequence)&eligible
            for k in keys:records[k].append(b[k][sel])
            #endfor
            n+=len(seq)
            if n%2000000<a.chunk:print(f'Particle entries scanned {n:,}; latest event_sequence {int(seq[-1]):,}',flush=True)
            #endif
            if np.any(~eligible):
                stopped_at_boundary=True
                print(f'Particle event boundary reached after scanning {n:,} particle rows.',flush=True)
                break
            #endif
        #endfor
        if not stopped_at_boundary and n>=a.max_particle_entries and n<t.num_entries:
            raise RuntimeError(f'Particle safety cap reached at {n:,} entries without reaching event_sequence {selected_last_sequence:,}. Increase --max-particle-entries only after checking event IDs.')
        #endif
    #endwith
    z={k:np.concatenate(v) for k,v in records.items()};proton=(z['gen_p']>0)&(z['rec_pindex']>=0)&(z['rec_pid']==2212)&(z['gen_t_proton']<0)&(z['rec_t_proton_raw']<0)
    fit=proton&(z['fitted_selected']>0)&(z['rec_t_proton_corrected']<0)
    print(f'Generated proton rows: {len(z["gen_p"]):,}; correctly PID-matched: {proton.sum():,}; fitted selected: {fit.sum():,}')
    fig,axs=plt.subplots(1,2,figsize=(13,5))
    bins=np.linspace(.10,.6,51)
    for sel,key,label in [(proton,'rec_t_proton_raw','Raw matched REC proton'),(fit,'rec_t_proton_corrected','Corrected fitted proton')]:
        x=-z['gen_t_proton'][sel];y=-z[key][sel]-x;idx=np.searchsorted(bins,x,side='right')-1;center=.5*(bins[:-1]+bins[1:]);mean=np.array([np.mean(y[idx==j]) if np.sum(idx==j)>=30 else np.nan for j in range(len(center))]);axs[0].plot(center,mean,'o-',ms=3,label=label)
    #endfor
    axs[0].axhline(0,color='black');axs[0].set(xlabel='Generated proton |t| (GeV²)',ylabel='Mean Δ|t| (GeV²)',title='MC::RecMatch: raw versus corrected');axs[0].legend();axs[0].grid(alpha=.2)
    v=proton&(z['gen_p']<1.5);axs[1].hist(z['rec_p'][v]-z['gen_p'][v],bins=120,range=(-.3,.3),histtype='step');axs[1].axvline(0,color='black');axs[1].set(xlabel='Matched raw proton Δp (GeV)',ylabel='Events',title='Correct PID and MC::RecMatch');fig.tight_layout();fig.savefig(out/'v8_03_truth_matched_protons.png',dpi=170);plt.close(fig)
    for lo,hi in [(.11,.15),(.15,.25),(.25,.4)]:
        sel=(g>=lo)&(g<hi);v=delta[sel];print(f'Event-level generated |t| {lo:.2f}-{hi:.2f}: N={len(v):,}, mean shift={np.mean(v):+.5f}, median={np.median(v):+.5f}, downward={np.mean(v<0):.1%}')
    #endfor
    print('WARNING: Candidate-level only. Event-level selection and MC::RecMatch proton sample differ; neither reproduces final 4D exclusivity cuts.')
    print('Output:',out)
if __name__=='__main__':main()
