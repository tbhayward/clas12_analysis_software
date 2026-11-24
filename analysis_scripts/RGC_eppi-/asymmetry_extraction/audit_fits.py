#!/usr/bin/env python3
"""Reconstruct all saved fits, audit positivity, score, curvature and starting values."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import nnls
from inputs import load,PERIODS
from likelihood import ConditionalLikelihood,NAMES
HERE=Path(__file__).resolve().parent

def audit(output,registry,metadata,dilution):
    groups,df,input_audit,identity=load(registry,metadata,dilution,HERE.parent/'channel_selection')
    original=json.loads((output/'input_audit.json').read_text())['identity']
    for key in ['cutset_id','metadata_sha256','dilution_sha256']:
        if original[key]!=identity[key]:raise ValueError('Audit input mismatch: '+key)
    folder=output/'audit';folder.mkdir(exist_ok=True)
    rows=[];details=[]
    for path in sorted((output/'fits').glob('*.json')):
        result=json.loads(path.read_text());b=result['bin'];scope=result['scope']
        names=PERIODS if scope=='combined' else [scope]
        selected=[]
        for name in names:
            g=groups[name];mask=g['events']['bin']==b
            selected.append(dict(dataset=name,events={k:v[mask] for k,v in g['events'].items() if k!='bin'},states=g['states']))
        if result['status']=='valid':
            beta=np.array(result['diluted_values']);saved_nll=result['nll']
        else:
            candidates=[c for c in result.get('candidates',[]) if c.get('values') is not None and np.isfinite(c.get('nll',np.nan))]
            if not candidates:
                rows.append(dict(scope=scope,bin=b,status=result['status'],primary_reason='missing_saved_parameters'));continue
            chosen=min(candidates,key=lambda c:c['nll']);beta=np.array(chosen['values']);saved_nll=chosen['nll']
        cost=ConditionalLikelihood(selected,target=len(beta)==7)
        numerator=1+cost.observed@beta;denominator=1+cost.average@beta;allowed=1+cost.allowed@beta
        gradient=np.sum(cost.average/denominator[:,None]-cost.observed/numerator[:,None],axis=0)
        oa=cost.observed/numerator[:,None];da=cost.average/denominator[:,None]
        hessian=oa.T@oa-da.T@da;eigen=np.linalg.eigvalsh(hessian)
        fisher_diag=np.sum(oa*oa,axis=0)
        score=float(np.max(np.abs(gradient)/np.sqrt(np.maximum(fisher_diag,1e-30))))
        active=allowed<1e-5
        kkt=None
        if active.any():
            try:
                multiplier,_=nnls(cost.allowed[active].T,gradient,maxiter=10000)
                kkt=float(np.linalg.norm(gradient-cost.allowed[active].T@multiplier)/max(np.linalg.norm(gradient),1e-30))
            except RuntimeError:pass
        worst=int(np.argmin(allowed));offset=0;provenance={}
        for g in selected:
            n=len(g['events']['phi'])
            for state in np.unique(g['states'][:,1:],axis=0):
                if offset<=worst<offset+n:
                    i=worst-offset;a=g['events']
                    provenance=dict(dataset=g['dataset'],event_index_in_bin=i,phi=float(a['phi'][i]),hPb=float(state[0]),Pt=float(state[1]),hPbPt=float(state[2]),observed_hPb=float(a['h'][i]*a['pb'][i]),observed_Pt=float(a['pt'][i]))
                    for key in ['runnum','evnum']:
                        if key in a:provenance[key]=int(a[key][i])
                offset+=n
        parts=cost.allowed[worst]*beta
        dominant=NAMES[int(np.argmin(parts))]
        covariance=result.get('diluted_covariance')
        minimum_sigma=None
        if covariance is not None:
            cv=np.array(covariance)
            variance=np.einsum('ij,jk,ik->i',cost.allowed,cv,cost.allowed)
            minimum_sigma=float(np.min(allowed/np.sqrt(np.maximum(variance,1e-30))))
        reasons=[]
        if allowed.min()<1e-5:reasons.append('positivity_boundary')
        if score>.1:reasons.append('nonstationary_unconstrained_score')
        if eigen[0]<=0:reasons.append('nonpositive_unconstrained_curvature')
        if minimum_sigma is not None and minimum_sigma<2:reasons.append('Gaussian_error_reaches_positivity_boundary')
        if result['status']!='valid':reasons.append('optimizer_or_covariance_invalid')
        if not reasons:reasons=['no_local_failure_found']
        nll=float(np.log(denominator).sum()-np.log(numerator).sum())
        if abs(nll-saved_nll)>1e-5:raise ValueError(f'{scope} bin {b}: saved likelihood not reproduced')
        all_nll=result.get('multistart_nll',[c['nll'] for c in result.get('candidates',[])])
        spread=float(max(all_nll)-min(all_nll)) if all_nll else None
        row=dict(scope=scope,bin=b,status=result['status'],events=cost.events,minimum_allowed_factor=float(allowed.min()),minimum_observed_factor=float(numerator.min()),active_constraints=int(active.sum()),score_standardized_max=score,hessian_min_eigenvalue=float(eigen[0]),hessian_max_eigenvalue=float(eigen[-1]),boundary_distance_in_gaussian_sigma=minimum_sigma,kkt_relative_residual=kkt,multistart_nll_spread=spread,dominant_cancelling_term=dominant,primary_reason='|'.join(reasons))
        rows.append(row)
        details.append(dict(**row,provenance=provenance,unit_plus_components={NAMES[i]:float(v) for i,v in enumerate(parts)},gradient=gradient.tolist(),hessian_eigenvalues=eigen.tolist(),beta=beta.tolist()))
    with (folder/'bin_audit.csv').open('w',newline='') as f:
        fields=list(dict.fromkeys(k for row in rows for k in row));w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)
    (folder/'bin_audit.json').write_text(json.dumps(dict(input_identity=identity,bins=details),indent=2,allow_nan=False)+'\n')
    nominal=[r for r in rows if r['scope']=='combined']
    text=['Nominal 24-bin audit','', 'Saved likelihoods reproduced from the supplied selected files and charge metadata.','']
    for row in nominal:
        text.append(f"Bin {row['bin']:2d}: {row['status']}; minimum allowed G={row.get('minimum_allowed_factor',float('nan')):.3g}; {row['primary_reason']}; dominant cancelling term={row.get('dominant_cancelling_term','unknown')}")
    # Controlled diagnostic: remove the two UU degrees of freedom, without
    # replacing the nominal result or borrowing amplitudes from another channel.
    from iminuit import Minuit
    fixed_checks=[]
    for b in range(1,25):
        selected=[]
        for name in PERIODS:
            g=groups[name];mask=g['events']['bin']==b
            selected.append(dict(events={k:v[mask] for k,v in g['events'].items() if k!='bin'},states=g['states']))
        cost=ConditionalLikelihood(selected,target=True)
        m=Minuit(cost,0.,0.,.02,.01,.01,.01,.01,name=NAMES,grad=cost.gradient)
        m.errordef=.5;m.fixed[0]=True;m.fixed[1]=True
        m.migrad(ncall=5000);m.hesse()
        check=dict(bin=b,diagnostic_only=True,UU_cos_phi=0,UU_cos_2phi=0,
                   valid=bool(m.valid and m.accurate),nll=float(m.fval),
                   minimum_state_factor=float(np.min(1+cost.allowed@np.array(m.values))))
        fixed_checks.append(check)
    (folder/'fixed_UU_zero_diagnostic.json').write_text(json.dumps(fixed_checks,indent=2,allow_nan=False)+'\n')
    text+=['','Controlled diagnostic with both UU cosine ratios fixed to zero (not a nominal replacement):',f"{sum(r['valid'] for r in fixed_checks)}/24 bins pass optimizer/covariance checks."]
    isolation=[]
    for b in range(1,25):
        for variant,fixed,excluded in [('UU_cos_phi_zero',[0],None),('UU_cos_2phi_zero',[1],None),('without_RGB_Fa19',[],PERIODS[1])]:
            selected=[]
            for name in PERIODS:
                if name==excluded:continue
                g=groups[name];mask=g['events']['bin']==b
                selected.append(dict(events={k:v[mask] for k,v in g['events'].items() if k!='bin'},states=g['states']))
            cost=ConditionalLikelihood(selected,target=True)
            attempts=[]
            for uu in [0.,.25,-.25]:
                start=[uu,-uu/2,.02,.01,.01,.01,.01]
                for j in fixed:start[j]=0.
                m=Minuit(cost,*start,name=NAMES,grad=cost.gradient);m.errordef=.5
                for j in fixed:m.fixed[j]=True
                m.migrad(ncall=4000);m.hesse();attempts.append(m)
            valid=[m for m in attempts if m.valid and m.accurate]
            best=min(valid or attempts,key=lambda m:m.fval)
            isolation.append(dict(bin=b,variant=variant,diagnostic_only=True,valid=bool(valid),minimum_state_factor=float(np.min(1+cost.allowed@np.array(best.values))),nll=float(best.fval)))
    (folder/'isolation_diagnostics.json').write_text(json.dumps(isolation,indent=2,allow_nan=False)+'\n')
    text+=['','Additional isolation checks (diagnostics only, not nominal fits):']
    for variant in ['UU_cos_phi_zero','UU_cos_2phi_zero','without_RGB_Fa19']:
        selected=[r for r in isolation if r['variant']==variant]
        text.append(f"{variant}: {sum(r['valid'] for r in selected)}/24 valid; remaining invalid bins: "+','.join(str(r['bin']) for r in selected if not r['valid']))
    text+=['','An active positivity constraint is not equivalent to an interior likelihood minimum. The large Minuit EDM near the hard invalid-likelihood penalty is not a reliable measure of physical error size. The direct analytic score and curvature above avoid that penalty. KKT residuals test whether the score can be balanced by nonnegative multipliers of the active positivity constraints. A small residual is evidence for a constrained boundary solution, not an unconstrained Gaussian fit.','', 'Boundary-distance diagnostics are local linear Gaussian checks, not confidence intervals. Reaching the boundary within two covariance standard deviations flags potentially misleading symmetric errors even for an optimizer-valid fit.','', 'A spin-conditional likelihood has no UU information when every spin-dependent amplitude is zero. Small spin signals can therefore give weak or nonregular UU estimates. This audit does not add external constraints or fix any amplitudes.']
    (folder/'audit_report.txt').write_text('\n'.join(text)+'\n')
    print('\n'.join(text[:8]))
    return rows

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=HERE/'output')
    p.add_argument('--selected-inputs',type=Path,default=HERE.parent/'channel_selection/output/selected_inputs.json')
    p.add_argument('--run-info',type=Path,default=HERE.parent/'imports/clas12_run_info_ND3_source.csv')
    p.add_argument('--dilution',type=Path,default=HERE.parent/'dilution_factor/output/joint_period_method1.json')
    a=p.parse_args();audit(a.output,a.selected_inputs,a.run_info,a.dilution)
