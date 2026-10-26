#!/usr/bin/env python3
"""Generated-proton angular coverage and event-level DVCS reconstruction response.

Uses paired MCEvents (one row per HIPO event). Does NOT equate geometric
acceptance with final DVCS selection, or a zero observed count with true zero
acceptance. No MC::Event.weight is used. Outputs PNGs and CSVs directly to --output.
"""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

T_EDGES = np.array([.11,.15,.25,.40,.60,.80,1.00])
THETA_EDGES = np.arange(0., 181., 2.)
P_EDGES = np.arange(0., 2.501, .025)


def ratio(num, den):
    return np.divide(num, den, out=np.full(np.shape(num), np.nan, dtype=float), where=den>0)


def write_csv(path, rows, fields):
    with open(path, 'w', newline='') as handle:
        w=csv.DictWriter(handle, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    #endwith


def load_cells(path):
    import re
    rows=[]
    with open(path) as handle:
        for line in list(handle)[2:]:
            x=re.split(r'\s+', line.strip())
            if len(x)>=12 and x[0].isdigit():
                rows.append((int(x[0]),float(x[4]),float(x[5]),float(x[7]),float(x[8]),float(x[10]),float(x[11])))
            #endif
        #endfor
    #endwith
    if not rows:
        raise RuntimeError('No 3D cells parsed from '+str(path))
    #endif
    return rows


def process_chunk(task):
    """Read one ROOT entry range in an independent process; return reduced counts."""
    input_path, start, stop, mandatory, cells, boundary = task
    import uproot
    with uproot.open(input_path) as root:
        b = root['MCEvents'].arrays(mandatory, entry_start=start, entry_stop=stop, library='np')
    good = ((b['gen_kinematics_valid'] > 0) & np.isfinite(b['gen_t1']) &
            (b['gen_t1'] < 0) & np.isfinite(b['gen_p1_p']) &
            (b['gen_p1_p'] > 0) & np.isfinite(b['gen_p1_theta']) &
            np.isfinite(b['gen_Q2']) & np.isfinite(b['gen_x']))
    shape1 = (len(T_EDGES)-1, len(THETA_EDGES)-1)
    shape2 = (len(P_EDGES)-1, len(THETA_EDGES)-1)
    empty1 = np.zeros(shape1, dtype=np.int64)
    empty2 = np.zeros(shape2, dtype=np.int64)
    counts = np.zeros((len(cells), 6), dtype=np.int64)
    if not np.any(good):
        return len(b['gen_t1']), 0, 0, empty1, empty1.copy(), empty2, empty2.copy(), counts
    g = -b['gen_t1'][good]
    p = b['gen_p1_p'][good]
    theta = np.degrees(b['gen_p1_theta'][good])
    q = b['gen_Q2'][good]
    x = b['gen_x'][good]
    rec = (b['rec_kinematics_valid'][good] > 0) & (b['rec_final_state'][good] > 0)
    tbin = np.searchsorted(T_EDGES, g, side='right')-1
    ib = np.searchsorted(THETA_EDGES, theta, side='right')-1
    pb = np.searchsorted(P_EDGES, p, side='right')-1
    valid = (tbin >= 0) & (tbin < shape1[0]) & (ib >= 0) & (ib < shape1[1])
    linear = tbin[valid]*shape1[1] + ib[valid]
    gen1 = np.bincount(linear, minlength=np.prod(shape1)).reshape(shape1)
    rec1 = np.bincount((tbin[valid & rec]*shape1[1]+ib[valid & rec]), minlength=np.prod(shape1)).reshape(shape1)
    valid2 = valid & (pb >= 0) & (pb < shape2[0])
    gen2 = np.bincount(pb[valid2]*shape2[1]+ib[valid2], minlength=np.prod(shape2)).reshape(shape2)
    rec2 = np.bincount(pb[valid2 & rec]*shape2[1]+ib[valid2 & rec], minlength=np.prod(shape2)).reshape(shape2)
    # Preselect by |t| before testing xB and Q2 cell boundaries.
    for k, (_, xmin, xmax, qmin, qmax, tmin, tmax) in enumerate(cells):
        v = (g >= tmin) & (g < tmax) & (x >= xmin) & (x < xmax) & (q >= qmin) & (q < qmax)
        if not np.any(v):
            continue
        counts[k] = (np.count_nonzero(v), np.count_nonzero(v & rec),
                     np.count_nonzero(v & (theta > boundary)),
                     np.count_nonzero(v & (theta < 40)),
                     np.count_nonzero(v & (theta >= 40) & (theta <= boundary)),
                     np.count_nonzero(v & (theta > boundary)))
    return len(b['gen_t1']), len(g), int(np.count_nonzero(rec)), gen1, rec1, gen2, rec2, counts


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', default='/volatile/clas12/thayward/dvcs_reprocessing_20261022/production/fa18_inb/dvcs_rec/50nA/all.mc_events.root')
    ap.add_argument('--output', default='/u/home/thayward/clas12_analysis_software/analysis_scripts/misc/diagnostic_dvcs')
    ap.add_argument('--binning-file', default='/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_bin_v2.csv')
    ap.add_argument('--max-events', type=int, default=1000000)
    ap.add_argument('--chunk', type=int, default=100000)
    ap.add_argument('--workers', type=int, default=4, help='Parallel ROOT readers; use 1 for serial execution')
    ap.add_argument('--theta-boundary', type=float, default=64.23, help='Illustrative polar-angle threshold, NOT a validated geometric cut')
    ap.add_argument('--min-generated', type=int, default=100, help='Minimum generated events for displayed 2D efficiency')
    args=ap.parse_args()
    if args.workers < 1 or args.chunk < 1: ap.error('--workers and --chunk must be positive')
    import uproot
    inp=Path(args.input)
    if inp.is_dir(): inp=inp/'all.mc_events.root'
    #endif
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    cells=load_cells(args.binning_file)
    nbin=len(T_EDGES)-1
    nt=len(THETA_EDGES)-1; npb=len(P_EDGES)-1
    generated=np.zeros((nbin,nt),dtype=np.int64)
    reconstructed=np.zeros_like(generated)
    generated_2d=np.zeros((npb,nt),dtype=np.int64)
    reconstructed_2d=np.zeros_like(generated_2d)
    # Cell columns: gen, reconstructed, theta above nominal, low-angle, middle-angle, high-angle
    cellcounts=np.zeros((len(cells),6),dtype=np.int64)
    n_processed=0; n_truth=0; n_reco=0
    with uproot.open(inp) as root:
        tree = root['MCEvents']
        mandatory = ['gen_kinematics_valid','gen_t1','gen_p1_p','gen_p1_theta','gen_Q2','gen_x','rec_kinematics_valid','rec_final_state']
        missing = [k for k in mandatory if k not in tree.keys()]
        if missing: raise RuntimeError('Missing required MCEvents branches: '+str(missing))
        total = min(args.max_events, tree.num_entries) if args.max_events > 0 else tree.num_entries
    tasks = ((str(inp), start, min(total, start+args.chunk), mandatory, cells, args.theta_boundary)
             for start in range(0, total, args.chunk))
    print(f'Reading {inp} entries {total:,} with {args.workers} worker(s)', flush=True)
    def accumulate(result):
        nonlocal n_processed, n_truth, n_reco
        processed, truth, reco, gen1, rec1, gen2, rec2, cc = result
        n_processed += processed; n_truth += truth; n_reco += reco
        generated[:] += gen1; reconstructed[:] += rec1
        generated_2d[:] += gen2; reconstructed_2d[:] += rec2
        cellcounts[:] += cc
        print(f'Processed {n_processed:,}/{total:,}', flush=True)
    if args.workers == 1:
        for task in tasks: accumulate(process_chunk(task))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            for result in executor.map(process_chunk, tasks):
                accumulate(result)
    # Summary by |t| bin
    rows=[]
    for k in range(nbin):
        n=int(generated[k].sum());nr=int(reconstructed[k].sum())
        # Angular boundary applied to binned theta gives approximate count; exact cell fraction separate.
        above=int(generated[k,THETA_EDGES[:-1]>=args.theta_boundary].sum())
        rows.append(dict(t_min=T_EDGES[k],t_max=T_EDGES[k+1],n_generated=n,n_reconstructed=nr,reconstruction_fraction=nr/n if n else np.nan,n_generated_theta_above_boundary_approx=above,fraction_above_boundary_approx=above/n if n else np.nan))
    #endfor
    write_csv(out/'dvcs_proton_theta_by_t_v1.csv',rows,list(rows[0]))
    cellrows=[]
    for k,c in enumerate(cells):
        idx,xmin,xmax,qmin,qmax,tmin,tmax=c
        n,nr,na,nfd,ncd,nout=map(int,cellcounts[k])
        cellrows.append(dict(cell=idx,xb_min=xmin,xb_max=xmax,q2_min=qmin,q2_max=qmax,t_min=tmin,t_max=tmax,n_generated=n,n_reconstructed=nr,reconstruction_fraction=nr/n if n else np.nan,n_theta_above_boundary=na,fraction_theta_above_boundary=na/n if n else np.nan,n_theta_below_40=nfd,n_theta_40_to_boundary=ncd,n_theta_above_boundary_repeat=nout))
    #endfor
    write_csv(out/'dvcs_proton_147_cells_v1.csv',cellrows,list(cellrows[0]))
    fig,ax=plt.subplots(figsize=(10,6))
    for k in range(nbin):
        y=generated[k].astype(float)
        if y.sum():ax.step(THETA_EDGES[:-1],y/y.sum(),where='post',label=f'{T_EDGES[k]:.2f}–{T_EDGES[k+1]:.2f} (GeV²)')
        #endif
    #endfor
    ax.axvline(40,color='black',ls=':',label='FD/CD split (40°)')
    ax.axvline(args.theta_boundary,color='black',ls='--',label=f'Illustrative boundary ({args.theta_boundary:g}°)')
    ax.set(xlabel='Generated proton polar angle (degrees)',ylabel='Fraction of generated events per 2° bin',title='Generated proton angular distribution by |t| bin')
    ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'dvcs_proton_01_generated_theta.png',dpi=175);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,6))
    for k in range(nbin):
        den=generated[k].astype(float);num=reconstructed[k].astype(float)
        y=ratio(num,den);y[den<args.min_generated]=np.nan
        ax.plot((THETA_EDGES[:-1]+THETA_EDGES[1:])/2,y,'o-',ms=2,label=f'{T_EDGES[k]:.2f}–{T_EDGES[k+1]:.2f}')
    #endfor
    ax.axvline(40,color='black',ls=':');ax.axvline(args.theta_boundary,color='black',ls='--')
    ax.set(xlabel='Generated proton polar angle (degrees)',ylabel='Fraction with reconstructed DVCS candidate',title=f'Event-level reconstruction fraction (≥{args.min_generated} generated per angular bin)',ylim=(0,1.05))
    ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'dvcs_proton_02_reconstruction_vs_theta.png',dpi=175);plt.close(fig)
    eff=ratio(reconstructed_2d,generated_2d);eff[generated_2d<args.min_generated]=np.nan
    fig,ax=plt.subplots(figsize=(11,6))
    m=ax.pcolormesh(THETA_EDGES,P_EDGES,eff,cmap='viridis',vmin=0,vmax=np.nanmax(eff) if np.any(np.isfinite(eff)) else 1,shading='auto')
    fig.colorbar(m,ax=ax,label='Reconstructed DVCS candidate / generated')
    ax.axvline(40,color='white',ls=':');ax.axvline(args.theta_boundary,color='white',ls='--')
    ax.set(xlabel='Generated proton polar angle (degrees)',ylabel='Generated proton momentum (GeV)',title='Candidate reconstruction fraction in analysis |t| range')
    ax.set_ylim(0,1.6);fig.tight_layout();fig.savefig(out/'dvcs_proton_03_p_theta_efficiency.png',dpi=175);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,6))
    cellx=np.arange(len(cellrows)); frac=np.array([r['fraction_theta_above_boundary'] for r in cellrows])
    ax.scatter(cellx,frac,s=12)
    ax.set(xlabel='Cell index in integrated_bin_v2.csv (ordered)',ylabel=f'Generated fraction θp > {args.theta_boundary:g}°',title='Generator angular extrapolation indicator across published 3D cells',ylim=(0,1.05))
    ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'dvcs_proton_04_cell_angular_fraction.png',dpi=175);plt.close(fig)
    print(f'Truth-valid events: {n_truth:,}; events with reconstructed DVCS candidate: {n_reco:,}')
    for row in rows:print(f"|t|={row['t_min']:.2f}–{row['t_max']:.2f}: generated={row['n_generated']:,}, candidate={row['n_reconstructed']:,}, above-boundary fraction≈{row['fraction_above_boundary_approx']:.2%}")
    #endfor
    print('IMPORTANT: reconstruction fraction is candidate-level, not final exclusive selection or detector-only proton tracking efficiency.')
    print('IMPORTANT: θ threshold is illustrative. A zero in a finite MC sample does NOT prove zero geometric acceptance.')
    print('IMPORTANT: these are generator-weighted fractions, not model-independent fractions of Nature.')
    print('Wrote files directly to',out)
#enddef

if __name__=='__main__':main()
