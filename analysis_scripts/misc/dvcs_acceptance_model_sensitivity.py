#!/usr/bin/env python3
"""Toy sensitivity to unknown cross-section shape beyond an illustrative proton angle.

Uses generated angular fractions from dvcs_proton_geometric_acceptance_v1.py.
Assumes zero reconstruction beyond the boundary, fixed accessible cross section,
and multiplies the inaccessible-region cross section by a factor r. This is a
sensitivity illustration, NOT a measured acceptance correction or systematic.
"""
import argparse,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',default='/u/home/thayward/clas12_analysis_software/analysis_scripts/misc/diagnostic_dvcs/dvcs_proton_147_cells_v1.csv')
    p.add_argument('--output',default='/u/home/thayward/clas12_analysis_software/analysis_scripts/misc/diagnostic_dvcs')
    p.add_argument('--relative-changes',default='-0.2,0.2',help='Fractional true inaccessible-region yield shifts relative to MC')
    p.add_argument('--min-generated',type=int,default=100)
    a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    with open(a.input,newline='') as f: rows=list(csv.DictReader(f))
    #endwith
    changes=[float(z) for z in a.relative_changes.split(',')]
    records=[]
    for row in rows:
        n=int(row['n_generated']);f=float(row['fraction_theta_above_boundary'])
        if n<a.min_generated or not np.isfinite(f):continue
        #endif
        r={'cell':int(row['cell']),'t_min':float(row['t_min']),'t_max':float(row['t_max']),'n_generated':n,'f_out':f}
        for d in changes:
            # If nature has r=(1+d) times the generated inaccessible yield,
            # sigma_extracted/sigma_true = Ngen_MC / Ngen_Nature.
            r[f'extracted_over_true_for_outside_shift_{d:+g}']=1/(1+d*f)
        #endfor
        records.append(r)
    #endfor
    if not records:raise RuntimeError('No cells with sufficient generated events')
    #endif
    with open(out/'dvcs_proton_model_sensitivity_v1.csv','w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    #endwith
    fig,ax=plt.subplots(figsize=(10,6))
    xx=np.array([r['cell'] for r in records])
    for d in changes:
        yy=np.array([100*(r[f'extracted_over_true_for_outside_shift_{d:+g}']-1) for r in records])
        ax.plot(xx,yy,'.',ms=4,label=f'Outside yield {d:+.0%} vs MC')
    #endfor
    ax.axhline(0,color='black',lw=.8)
    ax.set(xlabel='Published 3D cell number',ylabel='Illustrative extraction bias (%)',title='Toy model dependence from angularly inaccessible generated population')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'dvcs_proton_05_model_sensitivity.png',dpi=175);plt.close(fig)
    records.sort(key=lambda r:r['f_out'],reverse=True)
    print('Highest generated fractions beyond illustrative boundary:')
    for r in records[:15]:print(f"Cell {r['cell']:3d}, |t| {r['t_min']:.2f}-{r['t_max']:.2f}: f_out={r['f_out']:.1%}, Ngen={r['n_generated']:,}")
    #endfor
    print('Toy only: assumes zero acceptance outside boundary; not a measured systematic.')
    print('Output:',out)
#enddef
if __name__=='__main__':main()
