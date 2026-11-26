#!/usr/bin/env python3
"""Merge a validated release into the supplied Overleaf project, preserving other sections."""
import argparse,json,re,shutil,zipfile
from pathlib import Path


def figure(path,caption,label,width=r'\linewidth'):
 return '\n'.join([r'\begin{figure}[p]\centering',r'\includegraphics[width='+width+']{'+path+'}',r'\caption{'+caption+r'}\label{'+label+'}',r'\end{figure}',r'\clearpage'])+'\n'


def build(template,release,output,allow_preview=False):
 release=Path(release);output=Path(output);project=output/'overleaf';project.mkdir(parents=True,exist_ok=True)
 checklist=json.loads((release/'completion_checklist.json').read_text())
 if not checklist['period_checks_complete'] and not allow_preview:raise ValueError('Run the full release command with target-period fits before building the final note')
 with zipfile.ZipFile(template) as z:
  for name in z.namelist():
   p=(project/name).resolve()
   if project.resolve() not in p.parents and p!=project.resolve():raise ValueError('Unsafe ZIP member')
  z.extractall(project)
 shutil.copytree(release/'overleaf_additions',project,dirs_exist_ok=True)
 # Replace figures still referenced by existing selection/kinematic sections.
 for component,dest in [('kinematic_distributions','figures/channel_selection/kinematic_distributions'),('dilution_factor','figures/asymmetry_extraction/dilution_factor')]:
  for p in (release/'support'/component).rglob('*.png'):
   target=project/dest/p.relative_to(release/'support'/component/'figures') if 'figures' in p.relative_to(release/'support'/component).parts[:1] else project/dest/p.name
   target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
 groupmap={'RGB_spring2019':'rgb_sp19','RGB_fall2019':'rgb_fa19','RGB_spring2020':'rgb_sp20','RGC_summer2022':'rgc_su22','RGC_fall2022_spring2023':'rgc_fa22_sp23'}
 exclusive=project/'figures/channel_selection/exclusivity';exclusive.mkdir(parents=True,exist_ok=True)
 for group,prefix in groupmap.items():
  for p in (release/'support/channel_selection/fits'/group).glob('*.png'):shutil.copy2(p,exclusive/(prefix+'_'+p.name))
 # Find the exact historical group prefixes from figure references, supporting the original source naming.
 oldexclusive=(project/'sections/channel_selection/exclusivity.tex').read_text()
 for group in groupmap:
  alternatives={'RGB_spring2019':['rgb_sp19'],'RGB_fall2019':['rgb_fa19'],'RGB_spring2020':['rgb_sp20'],'RGC_summer2022':['rgc_su22','rgc_summer2022','su22'],'RGC_fall2022_spring2023':['rgc_fa22_sp23','rgc_fa22sp23','rgc_fasp','fa22_sp23']}[group]
  for prefix in alternatives:
   for p in (release/'support/channel_selection/fits'/group).glob('*.png'):
    if prefix+'_'+p.name in oldexclusive:shutil.copy2(p,exclusive/(prefix+'_'+p.name))
 meta=json.loads((release/'release_manifest.json').read_text());total=sum(r['rows'] for r in meta['datasets'].values())
 p=project/'main.tex';s=p.read_text();abstract=r'''This note presents preliminary beam-spin and longitudinal target-spin structure-function ratios for exclusive charged-pion electroproduction on the neutron, $en\to e'p\pi^-$, using unpolarized D$_2$ data from Run Group B and longitudinally polarized ND$_3$ data from Run Group C. The selected sample contains TOTAL candidate combinations. The beam-spin ratio is reported in 24 $(x_B,-t')$ bins, while each target-spin ratio is reported in 12 bins obtained by merging adjacent $-t'$ intervals. A period-conditional likelihood retains the run-dependent beam and target polarizations and charge normalization. Both unpolarized cosine ratios are fixed to zero, while all five polarized ratios float. The preliminary release reports $F_{LU}^{\sin\phi}/F_{UU}$, $F_{UL,\mathrm{lab}}^{\sin\phi}/F_{UU}$ and $F_{UL,\mathrm{lab}}^{\sin2\phi}/F_{UU}$ with statistical uncertainties. Comparisons with the ongoing Bhattacharya--Manley calculation are provided.'''.replace('TOTAL',f'{total:,}')
 s=re.sub(r'\\begin\{abstract\}.*?\\end\{abstract\}',lambda m:r'\begin{abstract}'+'\n'+abstract+'\n'+r'\end{abstract}',s,flags=re.S);p.write_text(s)
 binning=r'''\subsection{Binning}\label{sec:release_binning}
The four $x_B$ intervals are $[0.10,0.25)$, $[0.25,0.35)$, $[0.35,0.45)$ and $[0.45,0.60)$. The event selection and exclusivity calibration retain six $-t'$ intervals per $x_B$: $[0.05,0.25)$, $[0.25,0.45)$, $[0.45,0.65)$, $[0.65,0.85)$, $[0.85,1.05)$ and $[1.05,1.25)$ GeV$^2$. Bin numbers increase through $-t'$ first and then $x_B$.
The beam-spin result retains these 24 bins. For the two target-spin results, adjacent pairs of selected bins are merged into $[0.05,0.45)$, $[0.45,0.85)$ and $[0.85,1.25)$ GeV$^2$, yielding 12 bins. The target-sector likelihood is refitted to the merged events and the dilution factor is recalculated from summed control-target counts, rather than averaging fitted ratios or corrections. Original fine-bin exclusivity windows remain attached to each event. Points are displayed at the selected combined-sample mean $-t'$ in their respective bins.
'''
 (project/'sections/asymmetry_extraction/binning.tex').write_text(binning)
 # Describe only methods actually run for this release.
 method=r"""\subsection{Maximum Likelihood Estimation Method}\label{sec:MLE}
The extraction uses a conditional likelihood for the observed run and helicity state at the measured event kinematics. Each period has its own normalization. Let $Q_{rh}$ be the accepted integrated charge in run $r$ and helicity state $h$, and let $\mathcal G(\vec x;r,h)$ contain the polarized harmonics, run-dependent polarizations and event-dependent depolarization ratios. The conditional state probability is
\begin{equation}
 P(r,h\mid\vec x,p)=\frac{Q_{rh}\mathcal G(\vec x;r,h)}{\sum_{r',h'\in p}Q_{r'h'}\mathcal G(\vec x;r',h')}.
\end{equation}
The negative log likelihood sums $-\log P$ over all selected combinations in the bin. No binning in azimuth or artificial balancing of helicity samples is used. Detector acceptance cancels under the assumption that it is independent of the spin state at fixed kinematics within each period. Different periods may have different acceptances.
For RGC, measured helicity charges and signed run-dependent target polarization are retained. For RGB, the helicity scaler split is unavailable. The preliminary convention assigns $Q_{r+}=Q_{r-}=Q_r/2$ from the measured accepted total charge. This equal-exposure assumption is distinct from a measured helicity split. RGB has $P_t=0$. Source runs with no surviving exclusive candidate still contribute their accepted charge to the normalization.
Both unpolarized cosine ratios are fixed to zero. The five polarized ratios float simultaneously. The original 24-bin fit provides $LU\sin\phi$. Adjacent pairs of fine bins are merged and the event likelihood is refitted in twelve bins to provide the two $UL$ ratios. The two $LL$ amplitudes remain free nuisance amplitudes in both fits and are not shown in the conference results.
Target-dependent amplitudes are fitted initially as diluted products $\beta_j=D_f R_j$. A shared $D_f$ for the three RGC periods in each bin is obtained from the joint material fit. The covariance is transformed to $R_j$ and the first-order shared dilution contribution $\sigma_{D_f}^2(\partial R_i/\partial D_f)(\partial R_j/\partial D_f)$ is added. The statistical uncertainty is therefore conditional on the selected-sample and material model. It does not include unassigned selection or nuclear systematics.
Three starting values are tried. A fit is retained only with a valid minimum and accurate covariance. The cross-section factor must remain positive for every allowed spin state at the selected event kinematics. Saved fit files retain covariance matrices, correlations, optimizer diagnostics and positivity checks. Candidate combinations are treated as independent observations; repeated combinations from the same event are recorded in the sample audit.
A separate sensitivity study fixes $LL\cos\phi$ to the supplied BM stage-2 or stage-3 ratio, with both $UL$ ratios still free. Dilution propagation in that study includes the response of free amplitudes to the fixed physical $LL\cos\phi$ input. Neither constrained fit replaces the free-polarized baseline.
"""
 (project/'sections/asymmetry_extraction/mle.tex').write_text(method)
 p=project/'sections/asymmetry_extraction/pdf.tex';p.write_text(p.read_text().replace('All seven ratios are fitted simultaneously.','Both unpolarized cosine ratios are fixed to zero and all five polarized ratios are fitted simultaneously.'))
 (project/'sections/systematics/summary.tex').write_text(r"""\subsection{Preliminary Uncertainty Summary}
The displayed preliminary errors are statistical, including the shared dilution contribution for target-dependent ratios. Numerical point-to-point uncertainties from selection, nuclear backgrounds, bin migration and radiation have not yet been assigned. Beam and target polarization scale inputs are documented separately in the polarization section and are not included in the plotted error bars. The free-polarized extraction is retained as the baseline. Differences from the two BM-fixed $LL\cos\phi$ scenarios are sensitivity checks, not a complete systematic uncertainty estimate.
""")
 (project/'sections/systematics/run_period.tex').write_text(r"""\subsection{Run-period Comparisons}
Individual-period fits use the same nominal binning, frozen selection and shared dilution convention as the combined results. The comparisons are shown with the experimental results. Shared dilution errors correlate the three target-period results, and the periods populate different kinematics inside each bin. These comparisons are retained as diagnostics. No period-spread uncertainty is added to the preliminary error bars.
""")
 p=project/'sections/results/experimental.tex';s=r'\subsection{Experimental Results}\label{sec:release_results}'+'\n'+r'''The preliminary results use the free-polarized baseline described above. Only statistical errors are shown, including the propagated dilution uncertainty for the target-dependent ratios. Point-to-point selection, background, nuclear, radiative and polarization-scale uncertainties have not been assigned to these plots. The $LU$ result retains six $-t'$ bins per $x_B$ interval, while both $UL$ ratios use three. Individual points compatible with zero are retained. The $LL$ amplitudes float in the fits but are not included in the conference release.
\input{sections/generated/selected_counts}
'''
 s+=figure('figures/results/data_overlay_all_xB.png',r'Preliminary $LU$, $UL\sin\phi$ and $UL\sin2\phi$ ratios. Colors distinguish the four $x_B$ intervals. Error bars are statistical.','fig:release_overlay')
 for i in range(1,5):s+=figure(f'figures/results/data_xB_bin_{i}.png',rf'Preliminary spin-dependent ratios in $x_B$ interval {i}. The $LU$ panel has six $-t^\prime$ bins and each $UL$ panel has three.','fig:release_xb'+str(i))
 for obs in ['lu1','ul1','ul2']:s+='\n'+r'\input{sections/generated/'+obs+'_results_table}'+'\n'
 s+=r'''The period comparisons retain separate likelihood normalizations and the common dilution convention. They are diagnostic comparisons of samples with different within-bin kinematic populations. The two $UL$ observables are insensitive to the tested fixed-$LL\cos\phi$ scenarios relative to their statistical errors. The free baseline is therefore retained without a model constraint on $LL\cos\phi$.
'''
 for i in range(1,5):s+=figure(f'figures/results/periods_xB_bin_{i}.png',f'Period comparison in $x_B$ interval {i}. $LU$ includes the three RGB and three RGC samples; the target-spin panels show the three RGC ND$_3$ periods. Points are displaced horizontally for readability. Shared dilution uncertainties correlate target-period results.','fig:release_period'+str(i))
 p.write_text(s)
 p=project/'sections/results/mb_predictions.tex';s=r'''\subsection{Comparison to Bhattacharya--Manley Predictions}\label{sec:mb_comparison}
The preliminary Bhattacharya--Manley (BM) calculation includes longitudinal amplitudes with $\widetilde H_L$, $\widetilde E_L$ and the charged-pion pole, chiral-odd transverse contributions, and additional chiral-even transverse contributions. We show both the full BM prediction and cumulative stages: longitudinal/pole only, longitudinal/pole plus chiral-odd, and the full calculation. Each stage is a complete ratio with its own denominator. The curves are not additive components.
The supplied prediction tables use four fixed $x_B,Q^2$ settings and are not folded over the measured finite-bin kinematics. Their supplied lower and upper widths define the displayed bands and do not represent the full model uncertainty. Target-spin predictions are supplied with respect to the virtual-photon direction, whereas the extracted $UL$ ratios use the laboratory beam-axis target convention. They are overlaid directly in this preliminary comparison without a frame transformation. No GK curves are included for the neutron channel.
'''
 s+=figure('figures/results/BM_total_all_xB.png',r'Preliminary $\pi^-$ ratios compared with the full BM prediction. The $LU$ data retain six bins and the $UL$ data use three bins per $x_B$.','fig:release_bm_all')
 for i in range(1,5):s+=figure(f'figures/results/BM_stages_xB_bin_{i}.png',f'Cumulative BM comparisons in $x_B$ interval {i}. Orange, red and blue denote the longitudinal/pole, stage-2 and full predictions, respectively. Each curve is a separately normalized complete ratio.','fig:release_bm'+str(i))
 for obs in ['LU_sin_phi','UL_sin_phi','UL_sin_2phi']:s+=figure(f'figures/results/BM_stages_{obs}_2x2.png',r'Cumulative BM predictions and the preliminary data in all four $x_B$ intervals.','fig:release_bm_'+obs)
 p.write_text(s)
 # Update dilution numbers and state the released target binning explicitly.
 p=project/'sections/asymmetry_extraction/dilution.tex';s=p.read_text();consistency=json.loads((release/'support/dilution_factor/period_consistency.json').read_text())['summary'];xb,tp=consistency
 s=re.sub(r'Correlated statistical comparisons give .*?The comparisons include correlations from the shared Sp23 measurements\.',lambda m:f"Correlated statistical comparisons give $\\chi^2/\\mathrm{{ndf}}={xb['chi2']:.2f}/{xb['ndf']}$ ($p={xb['p_value']:.3f}$) for the $x_B$ dependence and ${tp['chi2']:.2f}/{tp['ndf']}$ ($p={tp['p_value']:.3f}$) for the $-t'$ dependence. The comparisons include correlations from the shared Sp23 measurements.",s,flags=re.S)
 start=s.find('For improved statistical stability')
 if start>=0:s=s[:start]+r'''For the released target-spin ratios, adjacent fine bins are merged before the material calculation. A shared-period dilution factor is fitted in each of the twelve merged bins, retaining independent period target rates and the common Sp23 CD$_2$/CH$_2$ transfer. Statistical errors are the standard deviations of 2,000 independent Poisson-replica fits per bin. The $LU$ ratio receives no direct dilution correction. The earlier 24-bin corrections remain available as a cross-check but are not used for the released three-bin $UL$ ratios.
\input{sections/generated/dilution_12_table}
'''
 s+=figure('figures/results/dilution_12_bins.png','Shared-period dilution factors in the twelve target-sector bins. Errors are the standard deviations from 2,000 Poisson-replica fits.','fig:release_dilution')
 p.write_text(s)
 # Remove stale detector-population numbers, preserving all distribution figures.
 p=project/'sections/channel_selection/kinematic_distributions.tex';s=p.read_text();start=s.find('Detection of the proton also');end=s.find(r'\FloatBarrier',start)
 if start>=0 and end>=0:s=s[:start]+r'''The detected-proton momentum and angular distributions are shown separately for FD and CD reconstruction. These are selected candidate populations rather than detector efficiencies. The current selected sample counts are given in Table~\ref{tab:release_counts}. The displayed distributions retain the fine selection bins, including for the merged target-spin analysis.
'''+s[end:]
 for sample in ['RGB','RGC']:
  for variable in ['Q2','W','phi']:
   s+=figure('figures/results/'+sample+'_'+variable+'_three_bin_2x2.png',sample+r' distributions in the three merged target-sector bins. Each curve is normalized over its displayed range. The legends give the $-t^\prime$ intervals in GeV$^2$.','fig:release_kin_'+sample+variable)
 p.write_text(s)
 p=project/'sections/channel_selection/pid.tex';p.write_text(r'''\subsection{Particle Identification}
The relaxed processing retains EventBuilder electron, proton and negative-pion candidates with raw electron momentum above 2~GeV. The offline selection requires finite necessary kinematics, all three fiducial selections, an FD pion, and an FD or CD proton, followed by the DIS and exclusivity cuts. The aggregate legacy processing cuts are not imposed in the offline selection. Consequently the relaxed sample is not identical to the earlier processing-level PID selection. The original scalar diagnostic fields remain available for selection and background studies.
''')
 # Old calibration tables contain obsolete fitted numbers. Keep refreshed figures, append current windows.
 p=project/'sections/channel_selection/exclusivity.tex';s=p.read_text();first=s.find(r'\begin{figure}')
 if first>=0:
  figures=s[first:]
  intro=r'''\subsection{Exclusivity Cuts}\label{sec:exclusivity}
The electron, proton and pion define the deuteron-reference missing mass $M_{X,d}^2=(k+P_d-k'-p_p-p_\pi)^2$ and full spatial missing momentum. Calibration is performed in the original 24 bins, separately for the three RGB periods and RGC Su22, with RGC Fa22+Sp23 pooled for window determination. Each RGC control target receives the same windows as ND$_3$ in its period. Period identities and luminosities remain separate.
A Gaussian spectator-proton peak plus nonnegative linear background defines the mass interval $\mu_M\pm2\sigma_M$. After this cut, a gamma signal plus higher-mode gamma tail describes missing momentum. The selected interval is the signal's central 95.45\% probability interval. A quadratic mass background, a varied tail shape and doubled histogram widths provide diagnostic comparisons. Unresolved Su22 shapes can be transferred from the Fa22+Sp23 reference with documented donor corrections. The signed pion-angle residual is a closure plot and is not an additional cut.
The relaxed sample is recalibrated before the three-bin target refit. All selected events retain their original fine-bin exclusivity windows. Current windows and any reference assistance are retained in the release metadata. The figures below show the new calibration, replacing the earlier sample's fits.
'''
  figures=''
  descriptions={'mass':r'deuteron-reference missing-mass-squared fits','momentum':r'missing-momentum fits after the mass cut','angle':r'pion polar-angle closure after both cuts'}
  for group,prefix in groupmap.items():
   for stage in ['mass','momentum','angle']:
    for xi in range(1,5):
     figures+=figure('figures/channel_selection/exclusivity/'+prefix+'_'+stage+'_xb'+str(xi)+'.png',group.replace('_',r'\_')+': '+descriptions[stage]+' in $x_B$ interval '+str(xi)+'. The panels retain the six fine selection bins.','fig:new_excl_'+prefix+'_'+stage+str(xi))
  p.write_text(intro+r'\input{sections/generated/exclusivity_windows}'+'\n'+figures)
 # Update the nominal RGB source run range.
 p=project/'sections/channel_selection/datasets.tex';s=p.read_text().replace('11093--11283','11093--11297');p.write_text(s)
 # Check all literal figure and input references before packaging.
 missing=[]
 for p in project.rglob('*.tex'):
  clean='\n'.join(line.split('%')[0] for line in p.read_text().splitlines())
  for name in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',clean):
   if not (project/name).exists():missing.append(name)
 if missing:raise ValueError('Missing referenced figures: '+', '.join(sorted(set(missing))))
 archive=output/('pi_minus_preliminary_release_overleaf.zip' if checklist['period_checks_complete'] else 'DRAFT_pipeline_preview_overleaf.zip')
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
  for p in project.rglob('*'):
   if p.is_file() and p.suffix not in ['.aux','.log','.out','.toc','.synctex','.pdf']:z.write(p,p.relative_to(project))
 print(archive)
 return project,archive


if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--template',type=Path,required=True);p.add_argument('--release',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--allow-preview',action='store_true',help='Allow a clearly named draft before period fits are available');a=p.parse_args();build(a.template,a.release,a.output,a.allow_preview)
