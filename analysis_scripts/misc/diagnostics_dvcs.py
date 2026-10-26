#!/usr/bin/env python3
"""CLAS12 paired-MC DVCS migration: global |t| and actual (xB,Q2,|t|) cells.
Unweighted MC::Event counting; independent even/odd response/closure; PNG/CSV only.
Extended from supplied diagnostics_dvcs(2).py. Candidate-level only (NOT full exclusivity/fiducial selection).
"""
import argparse, csv, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

EDGES=np.array([.11,.15,.25,.40,.60,.80,1.00])

def ratio(a,b):
    return np.divide(a,b,out=np.full(np.broadcast_shapes(np.shape(a),np.shape(b)),np.nan),where=np.asarray(b)!=0)

def unfold(y,R,prior,iters):
    t=np.maximum(prior.astype(float),1.e-12)
    eff=R.sum(axis=0)
    for _ in range(iters):
        pred=R@t
        t=np.where(eff>0,t*np.nan_to_num(R.T@np.nan_to_num(ratio(y,pred)))/np.maximum(eff,1.e-12),0)
    #endfor
    return t

def load_cells(path):
    cells=[]
    with open(path) as f:
        for line in f.readlines()[2:]:
            tok=[v.strip() for v in line.split('\t') if v.strip()]
            if len(tok)<11 or not tok[0].isdigit(): continue
            cells.append(dict(id=int(tok[0]),xb0=float(tok[2]),xb1=float(tok[3]),q0=float(tok[5]),q1=float(tok[6]),t0=float(tok[8]),t1=float(tok[9])))
        #endfor
    #endfor
    if not cells: raise ValueError('No cells parsed from binning CSV')
    return cells

def assign(q,x,t,cells):
    """Assign to exact published 3D cell; -1 denotes outside published bins."""
    ans=np.full(len(t),-1,dtype=np.int16)
    # Group by Q2/xB rectangular cells; only 147 small loops per chunk.
    groups={}
    for c in cells:
        key=(c['q0'],c['q1'],c['xb0'],c['xb1'])
        groups.setdefault(key,[]).append(c)
    #endfor
    for (q0,q1,x0,x1),group in groups.items():
        base=(q>=q0)&(q<q1)&(x>=x0)&(x<x1)
        if not np.any(base): continue
        for c in group:
            match=base&(t>=c['t0'])&(t<c['t1'])
            ans[match]=c['id']-1
        #endfor
    #endfor
    return ans

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',default='/volatile/clas12/thayward/dvcs_reprocessing_20261022/production/fa18_inb/dvcs_rec/50nA/all.mc_events.root')
    ap.add_argument('--output',default='/u/home/thayward/clas12_analysis_software/analysis_scripts/misc/diagnostic_dvcs')
    ap.add_argument('--binning-file',default='/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_bin_v2.csv')
    ap.add_argument('--max-events',type=int,default=10000000)
    ap.add_argument('--chunk',type=int,default=200000)
    ap.add_argument('--iterations',type=int,default=5)
    ap.add_argument('--slopes',default='0,-2,2')
    ap.add_argument('--min-cell-generated',type=int,default=200)
    args=ap.parse_args()
    import uproot
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    cells=load_cells(args.binning_file); C=len(cells)
    print(f'Loaded {C} published analysis cells from {args.binning_file}',flush=True)
    if sorted(c['id'] for c in cells)!=list(range(1,C+1)): raise ValueError('Nonsequential cell IDs')
    if not (abs(cells[0]['xb0']-0.062)<1e-5 and abs(cells[0]['q0']-1.0)<1e-5 and abs(cells[0]['t0']-0.11)<1e-5): raise ValueError('Binning CSV columns parsed incorrectly')
    n=len(EDGES)-1;m=n+2; slopes=[float(s) for s in args.slopes.split(',')]
    # Global |t| counts, split response/validation.
    gen=np.zeros((2,m)); mat=np.zeros((2,m,m)); pseudo_g=np.zeros((len(slopes),m));pseudo_r=np.zeros_like(pseudo_g)
    genall=np.zeros(m);matall=np.zeros((m,m)); origins=np.zeros((5,120));origin_edges=np.linspace(0,1.2,121)
    # Published-cell migration: rows reconstructed cell, columns generated cell.
    # Additional category C = outside all published analysis cells.
    G=C+1; cellgen=np.zeros((2,G)); cellmat=np.zeros((2,G,G));
    cellpseudo_g=np.zeros((len(slopes),G));cellpseudo_r=np.zeros_like(cellpseudo_g)
    cell_gen_all=np.zeros(G); cell_rec_all=np.zeros(G); cell_diag=np.zeros(G)
    cell_deltas=np.zeros(C);cell_delta2=np.zeros(C);cell_deltan=np.zeros(C)
    # Row-level origin categories: same cell, different t bin in same Q2/xB, different Q2/xB, outside.
    cell_orig=np.zeros((C,4))
    nraw=ntruth=nrec=0
    path=Path(args.input)
    if path.is_dir():
        files=sorted(path.glob('*.mc_events.root'))
        if len(files)!=1: raise ValueError(f'Expected one *.mc_events.root: {files}')
        path=files[0]
    #endif
    fields=['gen_kinematics_valid','rec_kinematics_valid','gen_t1','rec_t1','gen_Q2','rec_Q2','gen_x','rec_x','rec_final_state']
    with uproot.open(path) as f:
        tree=f['MCEvents']; missing=[v for v in fields if v not in tree.keys()]
        if missing: raise ValueError(f'Missing paired-MC branches: {missing}')
        total=min(args.max_events,tree.num_entries)
        for start in range(0,total,args.chunk):
            stop=min(total,start+args.chunk); a=tree.arrays(fields,entry_start=start,entry_stop=stop,library='np'); nraw+=stop-start
            valid=(a['gen_kinematics_valid']>.5)&np.isfinite(a['gen_t1'])&(a['gen_t1']<0)&np.isfinite(a['gen_Q2'])&np.isfinite(a['gen_x'])
            rec=(a['rec_kinematics_valid']>.5)&(a['rec_final_state']>.5)&np.isfinite(a['rec_t1'])&(a['rec_t1']<0)&np.isfinite(a['rec_Q2'])&np.isfinite(a['rec_x'])
            g=-a['gen_t1'][valid];r=-a['rec_t1'][valid];q=a['gen_Q2'][valid];x=a['gen_x'][valid];qr=a['rec_Q2'][valid];xr=a['rec_x'][valid];rr=rec[valid]
            gi=np.searchsorted(EDGES,g,side='right');ri=np.searchsorted(EDGES,r[rr],side='right')
            gc=assign(q,x,g,cells);rc=np.full(len(g),C,dtype=np.int16);tmp=assign(qr[rr],xr[rr],r[rr],cells);rc[rr]=np.where(tmp>=0,tmp,C)
            gc=np.where(gc>=0,gc,C)
            split=(np.arange(start,stop)[valid]%2).astype(int)
            ntruth+=len(g);nrec+=int(rr.sum())
            np.add.at(genall,gi,1);np.add.at(matall,(ri,gi[rr]),1)
            np.add.at(cell_gen_all,gc,1);np.add.at(cell_rec_all,rc[rr],1)
            for s in (0,1):
                keep=split==s; both=keep&rr
                np.add.at(gen[s],gi[keep],1);np.add.at(mat[s],(np.searchsorted(EDGES,r[both],side='right'),gi[both]),1)
                np.add.at(cellgen[s],gc[keep],1);np.add.at(cellmat[s],(rc[both],gc[both]),1)
            #endfor
            test=split==0
            for k,slope in enumerate(slopes):
                w=np.exp(np.clip(slope*(g[test]-EDGES[0]),-20,20));both=test&rr
                np.add.at(pseudo_g[k],gi[test],w)
                np.add.at(pseudo_r[k],np.searchsorted(EDGES,r[both],side='right'),np.exp(np.clip(slope*(g[both]-EDGES[0]),-20,20)))
                np.add.at(cellpseudo_g[k],gc[test],w)
                np.add.at(cellpseudo_r[k],rc[both],np.exp(np.clip(slope*(g[both]-EDGES[0]),-20,20)))
            #endfor
            for j in range(5):
                origins[j]+=np.histogram(g[rr&(np.searchsorted(EDGES,r,side='right')==j+1)],bins=origin_edges)[0]
            #endfor
            good=rr&(rc<C)
            recidx=rc[good];gidx=gc[good];d=r[good]-g[good]
            np.add.at(cell_deltas,recidx,d);np.add.at(cell_delta2,recidx,d*d);np.add.at(cell_deltan,recidx,1)
            # Same Q2/xB rectangle even when t bin differs, using boundaries of reconstructed cell.
            for c in cells:
                sel=good&(rc==c['id']-1)
                if not np.any(sel):continue
                same=(gc[sel]==c['id']-1)
                same_qx=(q[sel]>=c['q0'])&(q[sel]<c['q1'])&(x[sel]>=c['xb0'])&(x[sel]<c['xb1'])
                cat=np.where(same,0,np.where(same_qx,1,np.where(gc[sel]<C,2,3)))
                np.add.at(cell_orig,(np.full(cat.size,c['id']-1),cat),1)
            #endfor
            print(f'Processed {stop:,}/{total:,} events',flush=True)
        #endfor
    #endfor
    labels=['below']+[f'{EDGES[i]:.2f}–{EDGES[i+1]:.2f}' for i in range(n)]+['above']
    R=ratio(mat[1],gen[1][None,:]);R=np.nan_to_num(R)
    incoming=matall[1:n+1]; incoming_total=incoming.sum(axis=1)
    below=np.nan_to_num(ratio(incoming[:,0],incoming_total));same=np.nan_to_num(ratio(np.diag(incoming[:,1:n+1]),incoming_total));above=np.nan_to_num(ratio(incoming[:,-1],incoming_total));other=1-below-same-above
    fig,ax=plt.subplots(figsize=(9,7));img=np.ma.masked_where(matall==0,matall);im=ax.imshow(img,origin='lower',aspect='auto',norm=LogNorm(vmin=1,vmax=max(2,matall.max())));ax.set(xticks=np.arange(m),yticks=np.arange(m),xticklabels=labels,yticklabels=labels,xlabel='Generated |t| bin (GeV²)',ylabel='Reconstructed |t| bin (GeV²)',title='Global DVCS |t| migration');ax.tick_params(axis='x',rotation=40);fig.colorbar(im,ax=ax,label='Events');fig.tight_layout();fig.savefig(out/'01_t_migration_matrix.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));xx=np.arange(n);ax.plot(xx,ratio(matall.sum(axis=0)[1:-1],genall[1:-1]),'o-',label='Reconstructed anywhere');ax.plot(xx,ratio(np.diag(matall)[1:-1],genall[1:-1]),'s-',label='Same bin');ax.plot(xx,ratio(matall[1:-1,1:-1].sum(axis=0),genall[1:-1]),'^-',label='Inside analysis');ax.set(xticks=xx,xticklabels=labels[1:-1],xlabel='Generated |t| bin (GeV²)',ylabel='Probability',title='Truth-normalized reconstruction probabilities');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'02_reconstruction_efficiency.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));bottom=np.zeros(n)
    for vals,name in [(below,'Below analysis'),(same,'Same bin'),(other,'Other analysis bin'),(above,'Above analysis')]:ax.bar(xx,vals,bottom=bottom,label=name);bottom+=vals
    #endfor
    ax.set(xticks=xx,xticklabels=labels[1:-1],xlabel='Reconstructed |t| bin (GeV²)',ylabel='Fraction of reconstructed yield',ylim=(0,1.04),title='Origins of reconstructed events');ax.legend(fontsize=9);fig.tight_layout();fig.savefig(out/'03_incoming_migration_fractions.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));centers=(origin_edges[1:]+origin_edges[:-1])/2
    for j in range(5):
        if origins[j].sum():ax.step(centers,origins[j]/origins[j].sum(),where='mid',label=f'Rec {labels[j+1]}')
        #endif
    #endfor
    ax.axvline(EDGES[0],ls='--',c='black',label='Analysis threshold');ax.set(xlabel='Generated |t| (GeV²)',ylabel='Fraction per 0.01 GeV² bin',title='Generated origins of low-|t| reconstructed bins');ax.legend();fig.tight_layout();fig.savefig(out/'04_generated_origins_low_t.png',dpi=170);plt.close(fig)
    accept=np.nan_to_num(ratio(mat[1,1:n+1].sum(axis=1),gen[1,1:n+1]));closure=[]
    fig,axes=plt.subplots(len(slopes),1,figsize=(10,3.8*len(slopes)),squeeze=False)
    for k,slope in enumerate(slopes):
        truth=pseudo_g[k,1:n+1];conv=ratio(pseudo_r[k,1:n+1],accept);unf=unfold(pseudo_r[k],R,gen[1],args.iterations)[1:n+1]
        b1=100*(ratio(conv,truth)-1);b2=100*(ratio(unf,truth)-1)
        for j in range(n):closure.append(dict(slope=slope,bin=j+1,t_low=EDGES[j],t_high=EDGES[j+1],truth=truth[j],conventional_bias_pct=b1[j],matrix_bias_pct=b2[j]))
        #endfor
        ax=axes[k,0];ax.axhline(0,c='black',lw=1);ax.plot(xx,b1,'o-',label='Conventional');ax.plot(xx,b2,'s-',label=f'Bayesian ({args.iterations} iter)');ax.set(xticks=xx,xticklabels=labels[1:-1],ylabel='Closure bias (%)',title=f'Global closure: shape slope {slope:+g} (GeV⁻²)');ax.legend();ax.grid(alpha=.2)
        print(f'Global slope {slope:+g}: first-bin conventional {b1[0]:+.2f}%, matrix {b2[0]:+.2f}%')
    #endfor
    fig.tight_layout();fig.savefig(out/'05_generator_shape_closure.png',dpi=170);plt.close(fig)
    # Cell-level diagnostics: full reconstructed versus generated Q2/xB/t matching.
    cellR=np.nan_to_num(ratio(cellmat[1],cellgen[1][None,:]));cellaccept=np.nan_to_num(ratio(cellmat[1,:C,:].sum(axis=1),cellgen[1,:C]))
    cellclosure={}
    # High-dimensional response may be poorly populated: use response forward prediction for nominal and reweighted closure.
    # Iterative Bayes uses explicit outside category and zero-efficiency protection; DO NOT interpret sparse cells as validated.
    for k,slope in enumerate(slopes):
        truth=cellpseudo_g[k,:C];conv=ratio(cellpseudo_r[k,:C],cellaccept)
        unf=unfold(cellpseudo_r[k],cellR,cellgen[1],args.iterations)[:C]
        cellclosure[slope]=(100*(ratio(conv,truth)-1),100*(ratio(unf,truth)-1),truth)
        good=(cellgen[1,:C]>=args.min_cell_generated)&(truth>=args.min_cell_generated)&(cellmat[1,:C,:].sum(axis=1)>=30)
        print(f'Cell slope {slope:+g}: {good.sum()}/{C} cells meet minimum-count thresholds; median |conventional bias|={np.nanmedian(np.abs(cellclosure[slope][0][good])) if good.any() else np.nan:.2f}%')
    #endfor
    rows=[]
    for j,c in enumerate(cells):
        origins_total=cell_orig[j].sum();v=cell_orig[j];mean=cell_deltas[j]/cell_deltan[j] if cell_deltan[j] else np.nan
        rms=np.sqrt(max(0,cell_delta2[j]/cell_deltan[j]-mean**2)) if cell_deltan[j] else np.nan
        row=dict(cell=c['id'],xb_min=c['xb0'],xb_max=c['xb1'],Q2_min=c['q0'],Q2_max=c['q1'],t_min=c['t0'],t_max=c['t1'],generated=cell_gen_all[j],reconstructed=cell_rec_all[j],reconstructed_same_cell=cellmat[:,j,j].sum(),fraction_same_cell=v[0]/origins_total if origins_total else np.nan,fraction_other_t_same_Q2_xB=v[1]/origins_total if origins_total else np.nan,fraction_other_Q2_xB=v[2]/origins_total if origins_total else np.nan,fraction_outside_published_cells=v[3]/origins_total if origins_total else np.nan,delta_t_mean=mean,delta_t_rms=rms,nominal_response_generated=cellgen[1,j],nominal_response_reconstructed=cellmat[1,j,:].sum(),sufficient_statistics=int(cellgen[1,j]>=args.min_cell_generated and cellmat[1,j,:].sum()>=30))
        for slope in slopes:
            a,b,t=cellclosure[slope];tag=f'{slope:+g}'.replace('+','p').replace('-','m').replace('.','d')
            row[f'closure_conventional_{tag}_pct']=a[j];row[f'closure_matrix_{tag}_pct']=b[j]
        #endfor
        rows.append(row)
    #endfor
    def writecsv(name,records):
        with open(out/name,'w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
        #endfor
    #endfor
    writecsv('cell_migration_summary.csv',rows);writecsv('global_closure_summary.csv',closure)
    writecsv('global_bin_summary.csv',[dict(bin=j+1,t_low=EDGES[j],t_high=EDGES[j+1],generated=genall[j+1],reconstructed=incoming_total[j],fraction_below=below[j],fraction_same=same[j],fraction_other=other[j],fraction_above=above[j]) for j in range(n)])
    np.savetxt(out/'global_response_matrix.csv',R,delimiter=',');np.savetxt(out/'published_cell_response_matrix.csv',cellR,delimiter=',')
    def heatmap(vals,title,filename,vmin=None,vmax=None):
        fig,ax=plt.subplots(figsize=(13,5));im=ax.scatter([c['id'] for c in cells],[.5*(c['t0']+c['t1']) for c in cells],c=vals,s=35,cmap='viridis',vmin=vmin,vmax=vmax);fig.colorbar(im,ax=ax,label=title);ax.set(xlabel='Published (xB, Q², |t|) cell number',ylabel='Cell center |t| (GeV²)',title=title);fig.tight_layout();fig.savefig(out/filename,dpi=170);plt.close(fig)
    #endfor
    heatmap([r['fraction_same_cell'] for r in rows],'Fraction reconstructed from same published 3D cell','06_cell_retention_map.png',0,1)
    heatmap([r['fraction_other_Q2_xB'] for r in rows],'Fraction from different generated (Q²,xB) region','07_cross_Q2_xB_migration.png',0,1)
    heatmap([r['delta_t_mean'] for r in rows],'Mean |t|rec − |t|gen (GeV²)','08_t_resolution_bias.png')
    for slope in slopes:
        tag=f'{slope:+g}'.replace('+','p').replace('-','m').replace('.','d')
        heatmap(cellclosure[slope][0],'Conventional closure bias (%) slope '+str(slope),'09_cell_closure_conventional_'+tag+'.png',-30,30)
    #endfor
    print(f'Events read {nraw:,}; generated valid {ntruth:,}; reconstructed candidates {nrec:,}; published cells {C}')
    print(f'First global bin: {below[0]:.2%} below analysis, {other[0]:.2%} other analysis bins, {same[0]:.2%} same bin')
    print('WARNING: 3D cell unfolding can be prior-dependent and statistically unstable. Review per-cell counts and zero-slope closure before interpretation.')
    print('WARNING: candidate-level selection only; published-cell assignment uses actual generated AND reconstructed Q2, xB, t, but not phi or full exclusivity cuts.')
    print('MC::Event.weight is deliberately ignored.');print('PNG/CSV output:',out)

if __name__=='__main__':main()
