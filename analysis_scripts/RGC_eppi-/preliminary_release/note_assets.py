"""Portable Overleaf figures and generated numerical tables for the release."""
import json,shutil
from pathlib import Path


def make_assets(out,rows,means,meta,dfrows):
 note=out/'overleaf_additions';fig=note/'figures/results';fig.mkdir(parents=True,exist_ok=True)
 for p in (out/'figures').glob('*.png'):shutil.copy2(p,fig/p.name)
 s=note/'sections/generated';s.mkdir(parents=True,exist_ok=True)
 names={'lu1':r'$F_{LU}^{\sin\phi}/F_{UU}$','ul1':r'$F_{UL,\mathrm{lab}}^{\sin\phi}/F_{UU}$','ul2':r'$F_{UL,\mathrm{lab}}^{\sin2\phi}/F_{UU}$'}
 for obs in names:
  text=[r'\begin{longtable}{rrrrr}',r'\caption{'+names[obs]+r'. Statistical uncertainties include dilution propagation where applicable.}\label{tab:release_'+obs+r'}\\',r'\toprule',r'$x_B$ interval & $-t^\prime$ interval & $\langle-t^\prime\rangle$ & Ratio & Statistical error\\',r'\midrule\endfirsthead',r'\toprule $x_B$ interval & $-t^\prime$ interval & $\langle-t^\prime\rangle$ & Ratio & Statistical error\\\midrule\endhead']
  xb=[.10,.25,.35,.45,.60]
  for r in rows:
   if r['observable']!=obs:continue
   i=r['x_index'];text.append(f"{xb[i]:.2f}--{xb[i+1]:.2f} & {r['x_low']:.2f}--{r['x_high']:.2f} & {r['x']:.3f} & {r['value']:.4f} & {r['stat']:.4f}"+r'\\')
  text += [r'\bottomrule',r'\end{longtable}']
  (s/f'{obs}_results_table.tex').write_text('\n'.join(text)+'\n')
 counts=meta['datasets'];text=[r'\begin{table}[htbp]\centering',r'\begin{tabular}{lrr}\toprule',r'Sample & Selected combinations & Unique events\\\midrule']
 for name,r in counts.items():text.append(name.replace('_',r'\_')+f" & {r['rows']:,} & {r['unique_events']:,}"+r'\\')
 text += [r'\bottomrule\end{tabular}',r'\caption{Selected physics samples after the relaxed processing and offline selection.}\label{tab:release_counts}',r'\end{table}']
 (s/'selected_counts.tex').write_text('\n'.join(text)+'\n')
 text=[r'\begin{longtable}{rrrr}',r'\caption{Shared-period dilution factors in the twelve target-sector bins.}\\',r'\toprule Bin & $D_f$ & Statistical error & Joint deviance\\\midrule\endhead']
 for r in sorted(dfrows,key=lambda r:r['bin']):text.append(f"{r['bin']} & {r['df']:.4f} & {r['stat_std']:.4f} & {r['deviance']:.2f}"+r'\\')
 text += [r'\bottomrule\end{longtable}'];(s/'dilution_12_table.tex').write_text('\n'.join(text)+'\n')
 (out/'completion_checklist.json').write_text(json.dumps(dict(release_observables=['LU sinphi','UL sinphi','UL sin2phi'],LU_points=24,UL_points_each=12,LL_figures_in_release=False,period_checks_complete=meta['period_checks'],overleaf_next_step='Merge generated figures/tables and refresh narrative, cut summaries, dilution statistics, abstract, captions and BM discussion against the supplied source ZIP. Compile and visually inspect before packaging.'),indent=2)+'\n')

 cuts=json.loads((out/'support/channel_selection/nominal_cuts.json').read_text())
 text=[]
 for group,data in cuts['groups'].items():
  text += [r'\begin{longtable}{rrrrrl}',r'\caption{Current exclusivity windows for '+group.replace('_',r'\_')+r'. Mass-squared boundaries are in GeV$^2$ and missing-momentum boundaries are in GeV.}\\',r'\toprule Bin & $M^2_{\min}$ & $M^2_{\max}$ & $p_{\min}$ & $p_{\max}$ & Shape source\\\midrule\endhead']
  for b in data['bins']:
   if not b['active']:raise ValueError('Inactive selection bin in release')
   origin='reference-assisted' if b['reference_assisted'] else 'independent'
   text.append(f"{b['bin']} & {b['mass_low']:.4f} & {b['mass_high']:.4f} & {b['momentum_low']:.4f} & {b['momentum_high']:.4f} & {origin}"+r'\\')
  text += [r'\bottomrule\end{longtable}']
 (s/'exclusivity_windows.tex').write_text('\n'.join(text)+'\n')
