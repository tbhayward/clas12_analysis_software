"""Final editorial updates shared by the reproducible Overleaf builder."""
from pathlib import Path
import shutil,json,csv,re

def apply(project,release):
    r=Path(project)
    release=Path(release)
    for name,title,text in [
     ('bound_free_neutrons','Nuclear Effects and Residual Backgrounds',r'The stationary-neutron reconstruction does not remove Fermi motion, final-state interactions or residual nuclear and nonexclusive backgrounds. The detected proton and the missing-mass and missing-momentum selections provide constraints on the selected topology, but do not establish a background-free sample. No numerical correction or uncertainty from these effects is assigned in the preliminary results.'),
     ('bin_migration','Bin Migration',r'The preliminary results use reconstructed kinematics. A dedicated detector-response and nuclear-smearing study has not been used to assign migration uncertainties. No numerical bin-migration uncertainty is imported from the proton-channel analysis.'),
     ('radiative','Radiative Effects',r'No radiative correction or numerical radiative uncertainty is applied to this preliminary extraction. The corresponding proton-channel uncertainties are not substituted for a neutron-channel calculation.')]:
     (r/f'sections/systematics/{name}.tex').write_text(r'\subsection{'+title+'}\n'+text+'\n')
    p=r/'sections/introduction.tex';s=p.read_text().replace('preliminary seven-amplitude model','preliminary harmonic model').replace('the measured polarized structure-function ratios and laboratory-frame target- and beam--target-spin amplitudes','the measured beam-spin and laboratory-frame target-spin ratios');p.write_text(s)
    p=r/'sections/asymmetry_extraction/pdf.tex';s=p.read_text();start=s.index('The RGB target has');end=s.index('Both unpolarized cosine ratios are fixed',start)
    s=s[:start]+r'''The RGB target has $P_t=0$, so both target-dependent lines of Eq.~\ref{eq:cross_section_factor} vanish for RGB. The unpolarized cosine inputs are fixed to zero for all periods. The beam-spin ratio is shared between RGB and RGC within a fitted kinematic bin. The four laboratory-frame $UL$ and $LL$ ratios are constrained by the polarized RGC samples. Equality of the beam-spin ratio between the selected neutron samples remains a working assumption to be assessed against within-bin kinematic and background differences.
    
    '''+s[end:];p.write_text(s)
    p=r/'sections/asymmetry_extraction/depolarization.tex';s=p.read_text();s=s.replace('The reference study found little dependence on $-t$ and only a slight dependence on $x_B$. This is the working expectation here; the ND$_3$/RGB distributions will be checked before replacing event-level factors by any averages.','The extraction retains event-level factors throughout. No replacement by bin-averaged factors is made.');p.write_text(s)
    p=r/'sections/results/experimental.tex';s=p.read_text();at=s.index(r'\input{sections/generated/selected_counts}')
    text=r'''The beam-spin ratio shows a resolved change of sign with $-t'$ in the lowest $x_B$ interval. The value near $\langle-t'\rangle=0.54$~GeV$^2$ is $-0.1045\pm0.0213$, whereas the first bin is positive. At the highest $x_B$ the six values remain positive, approximately $0.08$--$0.10$. The intermediate intervals retain additional structure that would be lost by integrating over $-t'$.
    The three-bin $UL\sin\phi$ results resolve forward negative values in the two middle $x_B$ intervals: $-0.1174\pm0.0368$ and $-0.0912\pm0.0323$. Several higher-$-t'$ points are less precise. The $UL\sin2\phi$ measurements provide constraints with varying precision. The positive higher-$-t'$ values, including $0.4950\pm0.2285$ in $0.35<x_B<0.45$, are retained as measured points without treating them as established signals. Results compatible with zero remain part of the release.
    
    '''
    s=s[:at]+text+s[at:];s=s.replace('The two $UL$ observables are insensitive to the tested fixed-$LL\\cos\\phi$ scenarios relative to their statistical errors.','Across the tested BM-fixed scenarios, the largest changes in $UL\\sin\\phi$ and $UL\\sin2\\phi$ are 0.071 and 0.047 times the free-fit statistical error, respectively. The errors themselves are essentially unchanged. These are paired-fit sensitivity scales, not significances.');p.write_text(s)
    p=r/'sections/results/mb_predictions.tex';s=p.read_text();at=s.index(r'\begin{figure}')
    text=r'''The full prediction does not reproduce the negative higher-$-t'$ beam-spin values at low $x_B$, and it is above the positive data at the highest $x_B$. The target-spin comparisons have larger experimental uncertainties. Some forward negative $UL$ values are compatible with the curves, while positive higher-$-t'$ $UL\sin2\phi$ points differ in sign from the full prediction. The present precision and the stated kinematic and target-axis conventions limit the interpretation of those differences. Agreement or disagreement with a cumulative stage is not by itself an extraction of an individual GPD or GTMD moment.
    
    '''
    s=s[:at]+text+s[at:];p.write_text(s)
    # Current measured means are supplied both as a portable CSV and a LaTeX table.
    rows=list(csv.DictReader((release/'bin_kinematics.csv').open()))
    text=[r'\begin{longtable}{lrrrrr}',r'\caption{Measured combined-sample kinematic means in the reported bins. The two $UL$ ratios share the same twelve-bin means.}\\',r'\toprule Observable & Bin & $\langle x_B\rangle$ & $\langle Q^2\rangle$ & $\langle W\rangle$ & $\langle-t^\prime\rangle$\\\midrule\endhead']
    for obs,n,label in [('lu1',24,r'$LU$'),('ul1',12,r'$UL$')]:
     for b in range(1,n+1):
      d={x['variable']:float(x['mean']) for x in rows if x['observable']==obs and int(x['bin'])==b}
      text.append(label+f" & {b} & {d['x']:.3f} & {d['Q2']:.3f} & {d['W']:.3f} & {d['tp']:.3f}"+r'\\')
    text += [r'\bottomrule\end{longtable}']
    (r/'sections/generated/bin_kinematics_table.tex').write_text('\n'.join(text)+'\n')
    p=r/'sections/asymmetry_extraction/binning.tex';p.write_text(p.read_text()+'\n'+r'\input{sections/generated/bin_kinematics_table}'+'\n')
    data=r/'data';data.mkdir(exist_ok=True)
    for name in ['measurements.csv','bin_kinematics.csv','three_bin_histograms.csv','dilution_12_bins.json','period_results.csv','period_fit_summary.json','constraint_stability.json','release_manifest.json']:
     shutil.copy2(release/name,data/name)
    (r/'README.txt').write_text('Compile main.tex with pdfLaTeX, BibTeX, then pdfLaTeX twice.\nThis preliminary release reports fine-bin LU and three-bin UL sinphi and UL sin2phi. LL amplitudes float in the fits but have no release plots.\nThe data directory contains exact released values and provenance. Statistical errors include dilution where applicable. No point-to-point systematic errors have been assigned.\nSource: completed relaxed campaign release, 2026-11-26.\n')
    

    for name in ['binning','dilution']:
        p=r/f'sections/asymmetry_extraction/{name}.tex'
        text=p.read_text().replace('\\input{sections/generated/bin_kinematics_table}','\\clearpage\n\\input{sections/generated/bin_kinematics_table}').replace('\\input{sections/generated/dilution_12_table}','\\clearpage\n\\input{sections/generated/dilution_12_table}')
        p.write_text(text)
    p=r/'sections/results/experimental.tex'
    text=p.read_text()
    for obs in ['lu1','ul1','ul2']:
        text=text.replace('\\input{sections/generated/'+obs+'_results_table}','\\clearpage\n\\input{sections/generated/'+obs+'_results_table}')
    p.write_text(text)
    p=r/'sections/generated/exclusivity_windows.tex'
    p.write_text(p.read_text().replace('\\begin{longtable}','\\clearpage\n\\begin{longtable}'))
    shutil.copytree(release/'input_snapshots',r/'data/input_snapshots',dirs_exist_ok=True)
