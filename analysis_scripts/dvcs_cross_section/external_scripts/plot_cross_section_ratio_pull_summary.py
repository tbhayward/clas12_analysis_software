#!/usr/bin/env python3
"""Three-row cross-section ratio and signed-pull summaries from saved comparison tables.

No event reading or model evaluation. Pass-1 includes its separate normalization
reference band. Model panels do not assign that experimental allowance to theory.
"""
from pathlib import Path
import argparse,csv,json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
REGIONS=[('t','t_abs',r'$|t|$ (GeV$^2$)',[.11,.15,.25,.4,.6,.8,1.]),('xB','xB',r'$x_B$',[.062,.09,.118,.155,.204,.268,.357,.446,.581]),('Q2','Q2',r'$Q^2$ (GeV$^2$)',[1.,1.5,2.,3.,4.,6.,11.])]
def read(path):
    with path.open() as f:return list(csv.DictReader(f))
def source_points(directory,reference):
    if reference=='Pass-1':
        points=[]
        for x in read(directory/'tables/lee_hayward_exact_same_bin_points.csv'):
            points.append(dict(point_id=x['hayward_point_id'],xB=float(x['hayward_xB']),Q2=float(x['hayward_Q2']),t_abs=float(x['hayward_t_abs']),ratio=float(x['hayward_over_lee_direct']),pull=float(x['direct_pull_pointwise'])))
    else:
        points=[]
        for x in read(directory/'tables/pass2_model_publication_point_residuals.csv'):
            if x['model_name']!=reference or x['scenario']!='raw':continue
            data=float(x['data']);model=float(x['model']);unc=float(x['point_unc'])
            if not (data>0 and model>0 and unc>0):continue
            points.append(dict(point_id=x['point_id'],xB=float(x['xB']),Q2=float(x['Q2']),t_abs=float(x['t_abs']),ratio=data/model,pull=(data-model)/unc))
    if len({p['point_id'] for p in points})!=len(points):raise ValueError('Duplicate points in '+reference)
    if not points:raise ValueError('No supported points for '+reference)
    if not all(math.isfinite(p[k]) for p in points for k in ['xB','Q2','t_abs','ratio','pull']):raise ValueError('Nonfinite points for '+reference)
    return points

def summarize(points):
    result=[]
    for key,col,label,edges in REGIONS:
        for lo,hi in zip(edges[:-1],edges[1:]):
            p=[p for p in points if lo<=p[col]<hi]
            if not p:continue
            ratios=np.array([x['ratio'] for x in p]);pulls=np.array([x['pull'] for x in p])
            result.append(dict(variable=key,low=lo,high=hi,n=len(p),ratio_mean=np.mean(ratios),ratio_median=np.median(ratios),ratio_p16=np.percentile(ratios,16),ratio_p84=np.percentile(ratios,84),pull_mean=np.mean(pulls),pull_median=np.median(pulls),pull_rms=np.sqrt(np.mean(pulls**2))))
        if sum(x['n'] for x in result if x['variable']==key)!=len(points):raise ValueError('Projection excludes points in '+key)
    return result

def render(points,reference,out,normalization,dpi):
    summary=summarize(points)
    plt.rcParams.update({'font.size':11,'axes.labelsize':12,'xtick.labelsize':10,'ytick.labelsize':10})
    fig,axes=plt.subplots(3,2,figsize=(12,12))
    # A shared right-column scale in each figure, with space for large model pulls.
    pullmax=max(abs(x[k]) for x in summary for k in ['pull_mean','pull_median']);lim=max(2.,math.ceil(1.15*pullmax))
    for row,(key,col,label,edges) in zip(axes,REGIONS):
        values=[x for x in summary if x['variable']==key];x=np.arange(len(values));labels=[f"{v['low']:g}–{v['high']:g}" for v in values]
        for ax in row:
            ax.set_xticks(x,labels,rotation=30,ha='right');ax.set_xlabel(label);ax.grid(axis='y',alpha=.2)
            for spine in ax.spines.values():spine.set_visible(True)
        left,right=row
        if reference=='Pass-1':left.axhspan(1-normalization,1+normalization,color='#64a862',alpha=.18,label=f'Pass-1 ±{100*normalization:g}% normalization reference',zorder=0)
        med=np.array([v['ratio_median'] for v in values]);err=np.array([[v['ratio_median']-v['ratio_p16'] for v in values],[v['ratio_p84']-v['ratio_median'] for v in values]])
        left.errorbar(x,med,yerr=err,color='#245c9c',fmt='o-',capsize=3,lw=1.8,label='Median and central 68%',zorder=3)
        left.axhline(1,color='.25',lw=1);left.set_ylim(0,2.5);left.set_yticks(np.arange(0,2.51,.5));left.set_ylabel('Pass-2 / '+reference)
        # Mark percentile ranges clipped by the fixed display limit.
        for xx,v in zip(x,values):
            if v['ratio_p84']>2.5:left.plot(xx,2.47,marker='^',color='#245c9c',ms=5)
            if v['ratio_median']>2.5:left.annotate(f"median {v['ratio_median']:.2f}",(xx,2.47),xytext=(-5,-17),textcoords='offset points',ha='right',fontsize=8,color='#245c9c')
        right.plot(x,[v['pull_mean'] for v in values],'o-',color='#245c9c',lw=1.8,label='Mean pull')
        right.plot(x,[v['pull_median'] for v in values],'s--',color='#c62828',lw=1.8,label='Median pull')
        right.axhline(0,color='.25',lw=1);right.set_ylim(-lim,lim);right.set_ylabel('Signed pointwise pull')
    axes[0,0].legend(fontsize=8.5,loc='upper left');axes[0,1].legend(fontsize=9)
    title='Direct same-bin Pass-2 / Pass-1 comparison' if reference=='Pass-1' else 'Pass-2 / '+reference+' at the reporting kinematics'
    subtitle=f'Unshifted central values, {len(points):,} points, statistical and point-to-point uncertainties'
    fig.suptitle(title+'\n'+subtitle,fontsize=15,y=.98)
    fig.tight_layout(rect=[0,0,1,.94]);fig.subplots_adjust(hspace=.40)
    slug={'Pass-1':'pass1','BH':'bh','KM15':'km15'}[reference];path=out/(slug+'_ratio_pull_kinematic_summary.png');fig.savefig(path,dpi=dpi);plt.close(fig)
    with (out/(slug+'_ratio_pull_kinematic_summary.csv')).open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    return dict(reference=reference,points=len(points),figure=str(path),normalization_band_fraction=normalization if reference=='Pass-1' else None,right_axis_limit=lim)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--comparison-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,default=Path('output/ratio_pull_summaries'));p.add_argument('--pass1-normalization',type=float,default=.31);p.add_argument('--dpi',type=int,default=220);a=p.parse_args()
    if not 0<=a.pass1_normalization<1:p.error('normalization must be between zero and one')
    a.output_dir.mkdir(parents=True,exist_ok=True)
    manifest=[render(source_points(a.comparison_dir,r),r,a.output_dir,a.pass1_normalization,a.dpi) for r in ['Pass-1','BH','KM15']]
    (a.output_dir/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for x in manifest:print(x['reference'],x['points'],x['figure'])
if __name__=='__main__':main()
