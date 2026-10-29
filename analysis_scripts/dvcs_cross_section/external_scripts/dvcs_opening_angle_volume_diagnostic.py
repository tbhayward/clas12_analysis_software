#!/usr/bin/env python3
"""Standalone geometric and BH/KM15-weighted DVCS electron-photon opening-angle diagnostic.

Requires numpy, pandas, scipy, matplotlib. Optional: gepard for --model bh/km15.
Does not modify the input CSV. Phi is read in degrees and uses the same
phi -> pi - radians(phi) convention as the supplied Gepard integration.
"""
import argparse
import math
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.stats import qmc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MP = 0.938272
ME = 0.00051099895
COLS = ('xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max','phimin','phimax')

def four_vectors(x, q2, tabs, phi_deg, energy):
    """Return angle between outgoing e and gamma from exact 2->3 kinematics.

    Laboratory: incoming e along +z; scattered e in xz plane. Hadronic
    azimuth uses the Trento convention (phi is angle between lepton and
    hadron planes). Electron mass neglected consistently in DIS kinematics.
    """
    nu = q2/(2*MP*x)
    ep = energy-nu
    cth = 1-q2/(2*energy*ep)
    good = (ep>0)&(np.abs(cth)<=1)
    cth = np.clip(cth,-1,1)
    sth = np.sqrt(np.maximum(0,1-cth**2))
    # k'=(Ep, Ep*sin(theta), 0, Ep*cos(theta))
    qx = -ep*sth
    qz = energy-ep*cth
    qmag = np.hypot(qx,qz)
    q0 = nu
    # p'.q = (Q2 - t)/2, from (p+q-p')^2=0
    # proton final energy Ep'=Mp+nu-Egamma; from t=2Mp(Mp-Ep')
    pp0 = MP+tabs/(2*MP)
    pmag = np.sqrt(np.maximum(0,pp0**2-MP**2))
    # p'.q = Mp*nu + (t+Q2)/2 (since q'.p'=p.q + q^2? derive via q'=p+q-p')
    # q'.p= Mp*nu - tabs/2; E_gamma=nu-tabs/(2Mp)
    egamma = nu-tabs/(2*MP)
    # q.q'=(t-Q2)/2 = -(Q2+tabs)/2, and q'=q+p-p'
    # q dot q' = nu*Egamma - qvec dot gammavec = -(Q2+tabs)/2
    cosqg = (nu*egamma+(q2-tabs)/2)/(qmag*egamma)
    valid = good & (egamma>0) & (np.abs(cosqg)<=1+1e-8)
    cosqg = np.clip(cosqg,-1,1)
    sinqg = np.sqrt(np.maximum(0,1-cosqg**2))
    # q unit axis and perpendicular in electron scattering plane
    ux,uz = qx/qmag,qz/qmag
    # photon azimuth relative to lepton plane; sign irrelevant for e-gamma angle
    ph = np.deg2rad(phi_deg)
    gx = egamma*(cosqg*ux+sinqg*uz*np.cos(ph))
    gz = egamma*(cosqg*uz-sinqg*ux*np.cos(ph))
    costheta = (sth*gx+cth*gz)/egamma
    angle = np.rad2deg(np.arccos(np.clip(costheta,-1,1)))
    return angle,valid

def mask(x,q2,tabs,energy):
    y=q2/(2*MP*x*energy)
    w2=MP**2+q2*(1/x-1)
    eps2=4*MP**2*x*x/q2
    tmin=q2*(2*(1-x)*(1-np.sqrt(1+eps2))+eps2)/(4*x*(1-x)+eps2)
    return (y>0.19)&(y<0.8113919276)&(w2>4)&(tabs>=tmin)

def gepard_weights(x,q2,tabs,phi,energy,model):
    import gepard as g
    from gepard.fits import th_KM15
    weights=np.full(len(x),np.nan)
    for i in range(len(x)):
        try:
            pt=g.DataPoint(xB=float(x[i]),t=-float(tabs[i]),Q2=float(q2[i]),
                phi=math.pi-math.radians(float(phi[i])),observable='XS',
                frame='trento',process='ep2epgamma',exptype='fixed target',
                in1energy=float(energy),in1charge=-1,in1polarization=0,
                in2particle='p')
            pt.prepare()
            pref=float(th_KM15.PreFacSigma(pt))
            bh=pref*float(th_KM15.TBH2unp(pt))
            val=bh if model=='bh' else bh+pref*(float(th_KM15.TINTunp(pt))+float(th_KM15.TDVCS2unp(pt)))
            if np.isfinite(val) and val>=0: weights[i]=val
        except Exception:
            continue
    return weights

def evaluate(row,energy,n,angle_cut,seed,model):
    edges=np.array([float(row[c]) for c in COLS]).reshape(4,2)
    if np.any(edges[:,1]<=edges[:,0]): return None
    # Scrambled Sobol, same points for cut and uncut => low-variance ratio
    sampler=qmc.Sobol(d=4,scramble=True,seed=seed)
    points=qmc.scale(sampler.random_base2(int(math.log2(n))),edges[:,0],edges[:,1])
    x,q2,tabs,phi=points.T
    ok=mask(x,q2,tabs,energy)
    angle,valid=four_vectors(x,q2,tabs,phi,energy)
    ok &= valid
    cut=ok&(angle>angle_cut)
    cube=np.prod(edges[:,1]-edges[:,0])*math.pi/180 # phi in rad
    nominal=cube*ok.mean()
    selected=cube*cut.mean()
    result=dict(V_nominal=nominal,V_cut=selected,R_volume=selected/nominal if nominal>0 else np.nan,
                n_kin=int(ok.sum()),n_cut=int(cut.sum()),n_total=n)
    if model!='none':
        # Evaluate model on all kinematically valid points (no silent replacements)
        ids=np.flatnonzero(ok)
        weights=gepard_weights(x[ids],q2[ids],tabs[ids],phi[ids],energy,model)
        valid_model=np.isfinite(weights)
        result['n_model_valid']=int(valid_model.sum())
        result['n_model_invalid']=int((~valid_model).sum())
        if len(ids)>0 and valid_model.all() and np.sum(weights)>0:
            result['R_model']=float(np.sum(weights[angle[ids]>angle_cut])/np.sum(weights))
        else:
            result['R_model']=np.nan
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--csv',default='output/csvs/dvcs_pass2_analysis.csv')
    ap.add_argument('--output',default='output/opening_angle_diagnostic')
    ap.add_argument('--energies',nargs='+',type=float,default=[10.604,10.200])
    ap.add_argument('--angle',type=float,default=8.)
    ap.add_argument('--samples',type=int,default=4096,help='power of two per bin')
    ap.add_argument('--model',choices=['none','bh','km15'],default='none',help='Gepard model; can be slow')
    ap.add_argument('--max-rows',type=int,default=0)
    a=ap.parse_args()
    if a.samples<2 or a.samples&(a.samples-1): ap.error('--samples must be a power of two')
    if a.model!='none':
        try: import gepard
        except ImportError: ap.error('Install gepard to use --model bh/km15')
    df=pd.read_csv(a.csv,low_memory=False)
    missing=[c for c in COLS if c not in df]
    if missing: ap.error(f'Missing CSV columns: {missing}')
    out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    results=[]
    for energy in a.energies:
        count=0
        for idx,row in df.iterrows():
            try:
                phi_width=float(row['phimax'])-float(row['phimin'])
                if not 0<phi_width<359:continue
                vals=[float(row[c]) for c in COLS]
                if not np.all(np.isfinite(vals)):continue
            except (ValueError,TypeError):continue
            if a.max_rows and count>=a.max_rows:break
            r=evaluate(row,energy,a.samples,a.angle,13579+idx,a.model)
            if r is None:continue
            r.update(row_index=idx,beam_energy=energy,**{c:float(row[c]) for c in COLS})
            results.append(r);count+=1
            if count%100==0:print(f'[angle] E={energy:g}: {count} bins',flush=True)
    z=pd.DataFrame(results)
    if z.empty:raise SystemExit('No eligible phi-binned rows in CSV')
    z.to_csv(out/'bin_ratios.csv',index=False)
    with open(out/'summary.txt','w') as f:
        for energy,g in z.groupby('beam_energy'):
            v=g['R_volume'].dropna()
            lines=[f'Beam energy {energy:g} GeV: {len(g)} bins, {len(v)} valid',
                   f'Geometric ratio min={v.min():.6g}, median={v.median():.6g}, max={v.max():.6g}',
                   *[f'R_volume < {th}: {(v<th).sum()}' for th in (0.99,0.95,0.90)]]
            if 'R_model' in g:
                w=g['R_model'].dropna()
                lines += [f'{a.model} weighted: {len(w)} valid; min={w.min():.6g}, median={w.median():.6g}' if len(w) else 'No complete valid model bins']
            print('\n'.join(lines));f.write('\n'.join(lines)+'\n\n')
            fig,ax=plt.subplots(figsize=(7,4.5))
            ax.hist(v,bins=np.linspace(0,1,51))
            ax.set(xlabel=r'$V(\theta_{e\gamma}>8^\circ)/V_{\mathrm{nominal}}$',ylabel='Number of bins',title=f'{energy:g} GeV')
            fig.tight_layout();fig.savefig(out/f'volume_ratio_hist_{energy:g}.png',dpi=180);plt.close(fig)
            # Per 3D cell, phi dependence
            groups=g.groupby(['xBmin','xBmax','Q2min','Q2max','t_abs_min','t_abs_max'],dropna=False)
            for cell_id,(_,cell) in enumerate(groups):
                if len(cell)<2:continue
                cell=cell.sort_values('phimin')
                fig,ax=plt.subplots(figsize=(7,3.8))
                phi=(cell.phimin.to_numpy()+cell.phimax.to_numpy())/2
                ax.plot(phi,cell.R_volume.to_numpy(),'o-',label='Geometric')
                if 'R_model' in cell:ax.plot(phi,cell.R_model.to_numpy(),'s--',label=a.model.upper())
                ax.set(xlabel=r'$\phi$ (degrees)',ylabel='Surviving fraction',ylim=(-0.03,1.04),xlim=(0,360))
                ax.set_title(f'{energy:g} GeV; xB [{cell.xBmin.iloc[0]:.3f},{cell.xBmax.iloc[0]:.3f}], Q² [{cell.Q2min.iloc[0]:.2f},{cell.Q2max.iloc[0]:.2f}], |t| [{cell.t_abs_min.iloc[0]:.2f},{cell.t_abs_max.iloc[0]:.2f}]',fontsize=9)
                if 'R_model' in cell:ax.legend()
                fig.tight_layout();fig.savefig(out/f'phi_E{energy:g}_cell{cell_id:04d}.png',dpi=140);plt.close(fig)
    print(f'[angle] Outputs: {out.resolve()}')

if __name__=='__main__':main()
