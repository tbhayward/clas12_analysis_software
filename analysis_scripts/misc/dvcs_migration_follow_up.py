#!/usr/bin/env python3
"""Paired-MC DVCS follow-up: generated-conditioned t resolution, 3D flow balance,
neighbor origins, and correct unweighted acceptance-ratio statistical errors.
Candidate-level diagnostics only; does NOT reproduce final 4D DVCS selection.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

T_EDGES=np.array([.11,.15,.25,.40,.60,.80,1.00])

def load_cells(path):
    cells=[]
    with open(path) as f:
        for line in f.readlines()[2:]:
            v=[z.strip() for z in line.split('\t') if z.strip()]
            if len(v)>=10 and v[0].isdigit():
                cells.append(dict(id=int(v[0]),xb0=float(v[2]),xb1=float(v[3]),q0=float(v[5]),q1=float(v[6]),t0=float(v[8]),t1=float(v[9])))
            #endif
        #endfor
    #endfor
    if not cells or sorted(c['id'] for c in cells)!=list(range(1,len(cells)+1)):
        raise ValueError('Invalid or nonsequential analysis cells in binning CSV')
    #endif
    return cells

def assign(q,x,t,cells):
    out=np.full(len(t),-1,dtype=np.int16)
    for c in cells:
        mask=(q>=c['q0'])&(q<c['q1'])&(x>=c['xb0'])&(x<c['xb1'])&(t>=c['t0'])&(t<c['t1'])
        out[mask]=c['id']-1
    #endfor
    return out

def safe_ratio(a,b):
    return np.divide(a,b,out=np.full(np.broadcast_shapes(np.shape(a),np.shape(b)),np.nan),where=np.asarray(b)!=0)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',default='/volatile/clas12/thayward/dvcs_reprocessing_20261022/production/fa18_inb/dvcs_rec/50nA/all.mc_events.root')
    ap.add_argument('--binning-file',default='/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_bin_v2.csv')
    ap.add_argument('--output',default='/u/home/thayward/clas12_analysis_software/analysis_scripts/misc/diagnostic_dvcs')
    ap.add_argument('--max-events',type=int,default=10000000)
    ap.add_argument('--chunk',type=int,default=200000)
    ap.add_argument('--resolution-bins',type=int,default=90)
    ap.add_argument('--min-resolution-count',type=int,default=30)
    args=ap.parse_args()
    import uproot
    cells=load_cells(args.binning_file); C=len(cells); O=C
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    path=Path(args.input)
    if path.is_dir():
        found=sorted(path.glob('*.mc_events.root'))
        if len(found)!=1:raise ValueError(f'Expected one mc_events ROOT file, found {found}')
        path=found[0]
    #endif
    fields=['gen_kinematics_valid','rec_kinematics_valid','gen_t1','rec_t1','gen_Q2','rec_Q2','gen_x','rec_x','rec_final_state']
    # Matrix rows: reconstructed cells + outside, columns: generated cells + outside.
    # Missing reconstruction is NOT placed in the outside row; tracked separately.
    flow=np.zeros((C+1,C+1),dtype=np.int64)
    truth=np.zeros(C+1,dtype=np.int64)
    reconstructed=np.zeros(C+1,dtype=np.int64)
    rec_fail=np.zeros(C+1,dtype=np.int64)
    # Generated-conditioned fine-t resolution; generated region is inclusive of underflow.
    rt_edges=np.linspace(0,1.2,args.resolution_bins+1)
    rcount=np.zeros(args.resolution_bins,dtype=np.int64)
    rsum=np.zeros(args.resolution_bins);rsum2=np.zeros(args.resolution_bins)
    rpos=np.zeros(args.resolution_bins,dtype=np.int64);rneg=np.zeros_like(rpos)
    # Per 3D generated cell t resolution, plus per-reconstructed cell diagnostic.
    cgen_n=np.zeros(C,dtype=np.int64);cgen_d=np.zeros(C);cgen_d2=np.zeros(C)
    crec_n=np.zeros(C,dtype=np.int64);crec_d=np.zeros(C);crec_d2=np.zeros(C)
    total=valid_total=rec_total=0
    with uproot.open(path) as f:
        tree=f['MCEvents']
        missing=[x for x in fields if x not in tree.keys()]
        if missing:raise ValueError(f'Missing branches {missing}; available {list(tree.keys())}')
        stopall=min(args.max_events,tree.num_entries)
        for start in range(0,stopall,args.chunk):
            stop=min(start+args.chunk,stopall)
            a=tree.arrays(fields,entry_start=start,entry_stop=stop,library='np')
            valid=(a['gen_kinematics_valid']>.5)&np.isfinite(a['gen_t1'])&(a['gen_t1']<0)&np.isfinite(a['gen_Q2'])&np.isfinite(a['gen_x'])
            rec=(a['rec_kinematics_valid']>.5)&(a['rec_final_state']>.5)&np.isfinite(a['rec_t1'])&(a['rec_t1']<0)&np.isfinite(a['rec_Q2'])&np.isfinite(a['rec_x'])
            g=-a['gen_t1'][valid];r=-a['rec_t1'][valid]
            q=a['gen_Q2'][valid];x=a['gen_x'][valid]
            qr=a['rec_Q2'][valid];xr=a['rec_x'][valid]
            rr=rec[valid]
            gc=assign(q,x,g,cells);gc=np.where(gc>=0,gc,O)
            rc=np.full(len(g),O,dtype=np.int16)
            rc[rr]=np.where((tmp:=assign(qr[rr],xr[rr],r[rr],cells))>=0,tmp,O)
            np.add.at(truth,gc,1)
            np.add.at(rec_fail,gc[~rr],1)
            np.add.at(reconstructed,rc[rr],1)
            np.add.at(flow,(rc[rr],gc[rr]),1)
            d=r[rr]-g[rr]
            bi=np.searchsorted(rt_edges,g[rr],side='right')-1
            use=(bi>=0)&(bi<len(rcount))
            np.add.at(rcount,bi[use],1);np.add.at(rsum,bi[use],d[use]);np.add.at(rsum2,bi[use],d[use]**2)
            np.add.at(rpos,bi[use],(d[use]>0).astype(int));np.add.at(rneg,bi[use],(d[use]<0).astype(int))
            for idx,nv,sv,s2v in [(gc[rr],cgen_n,cgen_d,cgen_d2),(rc[rr],crec_n,crec_d,crec_d2)]:
                ok=idx<C
                np.add.at(nv,idx[ok],1);np.add.at(sv,idx[ok],d[ok]);np.add.at(s2v,idx[ok],d[ok]**2)
            #endfor
            total+=stop-start;valid_total+=len(g);rec_total+=int(rr.sum())
            print(f'Processed {stop:,}/{stopall:,} entries',flush=True)
        #endfor
    #endfor
    # Conventional acceptance A=N_rec_in_cell / N_gen_in_cell.
    # For Poisson generated events, covariance Cov(Nrec_i,Ngen_i)=N_ii.
    # Delta-method variance: Var(A)= (Nrec + A^2 Ngen - 2 A Nii)/Ngen^2.
    # For fixed total generated event count, this expression is also the leading multinomial ratio variance.
    rows=[]
    for i,c in enumerate(cells):
        ng=float(truth[i]); nr=float(reconstructed[i]); nii=float(flow[i,i])
        incoming=nr-nii;outgoing=ng-nii # includes nonreconstruction and out-of-analysis
        out_rec=float(flow[:,i].sum()-nii)
        out_missing=float(rec_fail[i])
        # These incoming/outgoing numbers are raw counts, NOT efficiency-corrected yields.
        A=nr/ng if ng else np.nan
        var=(nr+A*A*ng-2*A*nii)/(ng*ng) if ng else np.nan
        binomial_var=A*(1-A)/ng if ng and 0<=A<=1 else np.nan
        sameqx=0;diffqx=0;outside=flow[i,O]
        for j,other in enumerate(cells):
            if j==i:continue
            same=(other['q0']==c['q0'] and other['q1']==c['q1'] and other['xb0']==c['xb0'] and other['xb1']==c['xb1'])
            if same:sameqx+=flow[i,j]
            else:diffqx+=flow[i,j]
            #endif
        #endfor
        mn=cgen_d[i]/cgen_n[i] if cgen_n[i] else np.nan
        sd=np.sqrt(max(0,cgen_d2[i]/cgen_n[i]-mn*mn)) if cgen_n[i] else np.nan
        rows.append(dict(cell=c['id'],xb_min=c['xb0'],xb_max=c['xb1'],Q2_min=c['q0'],Q2_max=c['q1'],t_min=c['t0'],t_max=c['t1'],Ngen=int(ng),Nrec=int(nr),Nsame=int(nii),incoming=int(incoming),incoming_fraction=incoming/nr if nr else np.nan,outgoing_total=int(outgoing),outgoing_reconstructed_elsewhere=int(out_rec),outgoing_not_reconstructed=int(out_missing),net_incoming_minus_outgoing=int(incoming-outgoing),net_migration_reconstructed_only=int(incoming-out_rec),net_migration_over_Nrec=(incoming-out_rec)/nr if nr else np.nan,net_migration_over_Ngen=(incoming-out_rec)/ng if ng else np.nan,from_other_t_same_Q2_xB=int(sameqx),from_other_Q2_xB=int(diffqx),from_outside=int(outside),acceptance=A,acceptance_sigma=np.sqrt(max(0,var)) if np.isfinite(var) else np.nan,acceptance_rel_sigma=np.sqrt(max(0,var))/A if np.isfinite(var) and A>0 else np.nan,old_binomial_sigma=np.sqrt(binomial_var) if np.isfinite(binomial_var) else np.nan,generated_conditioned_delta_t_mean=mn,generated_conditioned_delta_t_std=sd,resolution_events=int(cgen_n[i])))
    #endfor
    with open(out/'v7_cell_flow_balance_and_acceptance_uncertainty.csv','w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    #endfor
    np.savetxt(out/'v7_response_counts.csv',flow,delimiter=',',fmt='%d',header='rows=reconstructed cells 1..C then outside; columns=generated cells 1..C then outside; nonreconstructed omitted')
    fine=[]
    for j in range(len(rcount)):
        mean=rsum[j]/rcount[j] if rcount[j] else np.nan
        std=np.sqrt(max(0,rsum2[j]/rcount[j]-mean**2)) if rcount[j] else np.nan
        fine.append(dict(t_gen_low=rt_edges[j],t_gen_high=rt_edges[j+1],n=int(rcount[j]),delta_t_mean=mean,delta_t_std=std,fraction_rec_t_below_gen=rneg[j]/rcount[j] if rcount[j] else np.nan))
    #endfor
    with open(out/'v7_generated_conditioned_t_resolution.csv','w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(fine[0]));writer.writeheader();writer.writerows(fine)
    #endfor
    fig,ax=plt.subplots(figsize=(10,5))
    xx=(rt_edges[:-1]+rt_edges[1:])/2
    enough=rcount>=args.min_resolution_count
    mu=safe_ratio(rsum,rcount);sigma=np.sqrt(np.maximum(0,safe_ratio(rsum2,rcount)-mu**2))
    ax.errorbar(xx[enough],mu[enough],yerr=sigma[enough]/np.sqrt(rcount[enough]),fmt='o',ms=3,capsize=2,label='Mean ± statistical SE')
    ax.plot(xx[enough],sigma[enough],'.-',label='Resolution standard deviation')
    ax.axhline(0,color='black',lw=1)
    for edge in T_EDGES:ax.axvline(edge,color='gray',alpha=.22,lw=.8)
    #endfor
    ax.set(xlabel='Generated |t| (GeV²)',ylabel='Reconstructed − generated |t| (GeV²)',title='Generated-conditioned proton-based |t| response')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'v7_01_generated_conditioned_t_resolution.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,5))
    ax.plot(xx[enough],safe_ratio(rneg,rcount)[enough],'o-')
    ax.axhline(.5,color='black',lw=1,ls='--')
    ax.set(xlabel='Generated |t| (GeV²)',ylabel='Fraction with |t|rec < |t|gen',ylim=(0,1),title='Direction of reconstructed |t| shifts, conditioned on generated |t|')
    ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'v7_02_t_shift_direction.png',dpi=180);plt.close(fig)
    # Low-|t| analysis cells in Q2/xB coordinates, with net migration and incoming fractions.
    for name,key,vmin,vmax in [('incoming_fraction','incoming_fraction',0,1),('net_reconstructed_migration','net_migration_over_Nrec',None,None),('acceptance_relative_uncertainty','acceptance_rel_sigma',0,None)]:
        fig,axs=plt.subplots(2,3,figsize=(16,9),sharex=True,sharey=True,constrained_layout=True)
        for j,ax in enumerate(axs.flat):
            selected=[z for z in rows if abs(z['t_min']-T_EDGES[j])<1e-5 and abs(z['t_max']-T_EDGES[j+1])<1e-5]
            if selected:
                vals=np.array([z[key] for z in selected],float)
                im=ax.scatter([.5*(z['xb_min']+z['xb_max']) for z in selected],[.5*(z['Q2_min']+z['Q2_max']) for z in selected],c=vals,s=75,cmap='viridis',vmin=vmin,vmax=vmax,edgecolors='black',linewidths=.2)
                fig.colorbar(im,ax=ax,shrink=.75)
            #endif
            ax.set_title(f'{T_EDGES[j]:.2f} < |t| < {T_EDGES[j+1]:.2f} (GeV²)')
            ax.set_xlabel('xB');ax.set_ylabel('Q² (GeV²)');ax.grid(alpha=.15)
        #endfor
        fig.suptitle(name.replace('_',' ').title())
        fig.savefig(out/f'v7_03_{name}_q2_xb.png',dpi=180);plt.close(fig)
    #endfor
    # Generated-bin outgoing versus reconstructed-bin incoming, explicitly split losses.
    low=[z for z in rows if abs(z['t_min']-.11)<1e-5 and abs(z['t_max']-.15)<1e-5]
    low=sorted(low,key=lambda z:z['Nrec'],reverse=True)[:12]
    fig,ax=plt.subplots(figsize=(12,5))
    xx=np.arange(len(low));ax.bar(xx-.2,[z['incoming'] for z in low],width=.4,label='Incoming reconstructed from other cells')
    ax.bar(xx+.2,[z['outgoing_reconstructed_elsewhere'] for z in low],width=.4,label='Generated here, reconstructed in other cells')
    ax.set(xticks=xx,xticklabels=[str(z['cell']) for z in low],xlabel='Published low-|t| cell ID',ylabel='Event count',title='Incoming versus outgoing reconstructed migration (highest-stat low-|t| cells)')
    ax.legend();fig.tight_layout();fig.savefig(out/'v7_04_low_t_incoming_outgoing.png',dpi=180);plt.close(fig)
    print(f'Raw events: {total:,}; generated valid: {valid_total:,}; reconstructed candidates: {rec_total:,}; published cells: {C}')
    print('Low-|t| cells (top reconstructed counts):')
    for z in low[:8]:
        print(f"  cell {z['cell']:3d}: Nrec={z['Nrec']:6d}, incoming={z['incoming_fraction']:.1%}, net rec migration={z['net_migration_reconstructed_only']:+d}, A={z['acceptance']:.4g}, sigma_A={z['acceptance_sigma']:.3g}")
    #endfor
    print('IMPORTANT: Net reconstructed migration = incoming - reconstructed elsewhere; failed reconstruction is reported separately.')
    print('IMPORTANT: Statistical uncertainty uses covariance Cov(Nrec_i,Ngen_i)=Nsame_i, not the invalid subset-binomial expression.')
    print('LIMITATION: Candidate-level response only, not final 4D acceptance, exclusivity, phi, or detector-fiducial selections.')
    print('MC::Event.weight ignored. Generated particles have no detector topology; no FD/CD assignment made.')
    print('Output:',out.resolve())

if __name__=='__main__':main()
