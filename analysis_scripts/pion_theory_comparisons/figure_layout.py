"""Human-readable figure names and a compact, pi+-first gallery."""
from pathlib import Path
import html
OBS_NAMES={'u1':'UU_cos_phi','u2':'UU_cos_2phi','lu1':'LU_sin_phi','ul1':'UL_sin_phi','ul2':'UL_sin_2phi','ll0':'LL','ll1':'LL_cos_phi'}
OBS_TEXT={'u1':'UU cosφ','u2':'UU cos2φ','lu1':'LU sinφ','ul1':'UL sinφ','ul2':'UL sin2φ','ll0':'LL','ll1':'LL cosφ'}

def figure_path(output,old_path):
 rel=Path(old_path).relative_to(output);parts=rel.parts
 if parts[0]=='focused':
  suffix={'LU_pip_pim':'pi_plus_pi_minus_LU_sin_phi_2x2','LU_pip_pim_with_pi0':'pi_plus_pi_minus_LU_sin_phi_with_pi0_theory_2x2','LU_pip_pim_with_pi0_combined_errors':'pi_plus_pi_minus_LU_sin_phi_with_pi0_theory_combined_errors_2x2','LU_pip_pim_residuals':'pi_plus_pi_minus_LU_sin_phi_pointwise_residuals_2x2'}.get(rel.stem,rel.stem)
  return Path(output)/'figures/charge_comparison'/(suffix+'.png')
 species,variant=parts[:2];channel={'pip':'pi_plus','pim':'pi_minus'}[species]
 base=Path(output)/'figures'/channel
 if variant=='combined_errors':base=base/'alternative_error_display'
 kind=parts[2];stem=rel.stem
 if kind=='focused':
  obs=next(k for k in OBS_NAMES if stem.startswith(k+'_'));mode=stem[len(obs)+1:]
  obsname=OBS_NAMES[obs]
  if species=='pip' and obs in ['ul1','ul2'] and mode.startswith('cff_') and variant=='stat_syst':
   folder=base/'featured_stories'
  elif mode.startswith('cff_'):folder=base/'cff_breakdowns'/('all_products' if mode=='cff_all' else 'leading_products')/'four_xB_panels'
  else:folder=base/'total_comparisons/four_xB_panels'
  description={'cff_all':'all_CFF_products','cff_leading5':'leading_CFF_products','totals':'BM_total'}[mode]
  return folder/f'{channel}_{obsname}_{description}_2x2.png'
 xb='all_xB' if stem=='all_xB' else 'xB_'+stem.split('_')[-1].zfill(2)
 if kind=='totals':folder=base/'total_comparisons/seven_ratio_canvases';desc='BM_totals'
 elif kind=='bm_and_gk':folder=base/'BM_and_GK/seven_ratio_canvases';desc='BM_and_GK'
 else:folder=base/'cff_breakdowns'/('all_products' if kind=='cff_all' else 'leading_products')/'seven_ratio_canvases';desc='all_CFF_products' if kind=='cff_all' else 'leading_CFF_products'
 return folder/f'{channel}_seven_ratios_{desc}_{xb}.png'

def legacy_path(output,old_relative):
 p=Path(old_relative);parts=p.parts
 variant='alternative_error_display' if 'combined_errors' in parts else 'statistical_and_systematic'
 kind='BM_and_GK' if 'both_theories' in parts else 'BM_only' if 'bhattacharya_manley' in parts else 'GK_only'
 suffix='all_xB' if 'all_xB' in p.stem else 'xB_'+p.stem.split('_')[-1].zfill(2)
 return Path(output)/'figures/pi_plus/reference_original_layouts'/variant/kind/f'pi_plus_original_layout_{kind}_{suffix}.png'

def friendly_title(path):
 stem=Path(path).stem
 if stem.startswith('pi_plus_pi_minus'):return 'π+ and π− — '+('UL sin2φ' if 'UL_sin_2phi' in stem else 'UL sinφ' if 'UL_sin_phi' in stem else 'LU sinφ')+(' (pointwise residuals)' if 'residuals' in stem else ' (with π0 theory)' if 'with_pi0' in stem else '')
 for key,text in sorted(OBS_NAMES.items(),key=lambda kv:-len(kv[1])):
  if '_'+text+'_' in stem:return ('π+ ' if stem.startswith('pi_plus_') else 'π− ')+OBS_TEXT[key]+(' — summed chiral sectors' if 'chiral_sector_sums' in stem else ' — individual products colored by sector' if 'products_by_chiral_sector' in stem else ' — every supplied CFF product' if 'all_CFF' in stem else ' — leading products + remainder' if 'leading_CFF' in stem else ' — BM total')
 return stem.replace('pi_plus_pi_minus','π+ and π−').replace('pi_plus','π+').replace('pi_minus','π−').replace('_',' ')

def write_gallery(output,entries):
 output=Path(output)
 paths=sorted({e['path'] for e in entries})
 def cards(items):
  return '<div class="grid">'+''.join('<article><a href="'+html.escape(p)+'"><img loading="lazy" src="'+html.escape(p)+'"></a><h3>'+html.escape(friendly_title(p))+'</h3><small>'+html.escape(p)+'</small></article>' for p in items)+'</div>'
 featured=[p for p in paths if '/step_predictions/bm/' in p and '_2x2.png' in p and p.startswith('figures/pi_plus/')] or [p for p in paths if '/featured_stories/' in p]
 charge=[p for p in paths if '/charge_comparison/' in p]
 body='''<!doctype html><html><head><meta charset="utf-8"><title>Pion theory comparisons</title><style>body{font:16px system-ui;color:#202b33;background:#f4f5f7;max-width:1500px;margin:35px auto;padding:0 25px}h1{font-size:36px;margin-bottom:8px}p{line-height:1.6}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:22px}article{background:white;padding:16px;border:1px solid #dce1e5;border-radius:8px}img{width:100%;height:390px;object-fit:contain}h3{font-size:17px;margin:8px 0}small{font-size:11px;color:#667681;overflow-wrap:anywhere}h2{margin-top:42px}details{background:white;border:1px solid #ddd;padding:14px;margin:12px 0}summary{cursor:pointer;font-weight:600}.note{background:#e9f0f5;padding:15px;border-radius:7px}a{color:#176397}nav a{margin-right:24px}@media(max-width:850px){.grid{grid-template-columns:1fr}}</style></head><body><h1>Pion theory comparisons</h1><p>Updated BM totals and CFF products, with π+ first. Data points are unchanged. Curves use the authors’ supplied central values.</p><nav><a href="#featured">π+ cumulative predictions</a><a href="#plus">All π+ figures</a><a href="#minus">π− figures</a><a href="#charge">Charge comparison</a></nav><div class="note"><b>How to read the cumulative predictions:</b> orange is longitudinal/pole only, red adds chiral-odd contributions, and blue is the full BM model including chiral-even transverse contributions. Each is a complete ratio with its own denominator. The curves are not additive. Light bands show supplied theory fit uncertainties. Light-gray curves are GK where available. The separate CFF-product plots retain their original common-denominator interpretation. π+ has statistical bars and point-to-point systematic boxes. π− LU uses six bins and UL/LL use two. Its fixed UU inputs are not measurements.</div>'''
 sector_sums=[p for p in featured if 'chiral_sector_sums' in p]
 sector_sums.sort(key=lambda p:next((i for i,o in enumerate(['UL_sin_phi','UL_sin_2phi','LU_sin_phi','LL','LL_cos_phi','UU_cos_phi','UU_cos_2phi']) if '_'+o+'_' in p),99))
 body+='<h2>Shared plots — seven π+ ratios and π− LU sinφ</h2>'
 for variant in ['bm','bm+gk','bm+gk_combined_errors']:
  body+='<details><summary>'+('BM sectors' if variant=='bm' else 'BM + GK, combined statistical and point-to-point errors' if variant=='bm+gk_combined_errors' else 'BM sectors + light-gray GK curves')+'</summary><p><a href="figures/shared_plots/'+variant+'/pion_chiral_sector_comparisons_9_pages.pdf">Nine-page PDF</a></p>'+cards([p for p in paths if p.startswith('figures/shared_plots/'+variant+'/')])+'</details>'
 body+='<h2 id="featured">π+ cumulative predictions — all seven ratios</h2>'+cards(sector_sums or featured)
 if sector_sums:
  details=[p for p in featured if p not in sector_sums]
  body+='<details><summary>Inspect individual CFF products and the earlier detailed views</summary>'+cards(details)+'</details>'
 for species,title,anchor in [('pi_plus','π+ complete figure collection','plus'),('pi_minus','π− complete figure collection','minus')]:
  body+='<h2 id="'+anchor+'">'+title+'</h2>'
  for folder,label in [('step_predictions','Cumulative BM models, each with its own denominator'),('total_comparisons','Full BM totals with data'),('cff_breakdowns/chiral_sectors','CFF breakdowns — longitudinal/pole, chiral-odd and added chiral-even sectors'),('cff_breakdowns/leading_products','CFF breakdowns — leading products + remainder'),('cff_breakdowns/all_products','CFF breakdowns — every supplied product'),('BM_and_GK','BM and GK totals'),('alternative_error_display','Alternative π+ error display: statistical + systematic in quadrature'),('reference_original_layouts','Reference: the original π+ plotting layouts')]:
   items=[p for p in paths if p.startswith('figures/'+species+'/'+folder+'/')]
   if items:body+='<details><summary>'+label+' ('+str(len(items))+')</summary>'+cards(items)+'</details>'
 body+='<h2 id="charge">π+ / π− LU comparison</h2>'+cards(charge)
 sensitivity=[p for p in paths if p.startswith('figures/cff_sensitivity/')]
 if sensitivity:body+='<h2>CFF contribution-size study</h2>'+cards(sensitivity)
 body+='<p>Theory has not been folded over the experimental bins or transformed to the lab target axis. No updated theory error bands were supplied. See README.md for conventions and regeneration commands.</p></body></html>'
 (output/'index.html').write_text(body)
