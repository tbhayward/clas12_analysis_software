"""Bold Brandon-inspired colors, chosen consistently within each observable.

The new input is individual normalized CFF products. We do not relabel these as
Brandon's earlier grouped A/P curves or reconstruct pieces hidden in 'other'.
"""
import numpy as np
BLUE='#1f77b4';ORANGE='#ff7f0e';GREEN='#2ca02c';RED='#d62728';PURPLE='#9467bd'
COLORS=[BLUE,ORANGE,GREEN,RED,PURPLE,'#17becf','#8c564b','#e377c2','#7f7f00','#165e70','#a55194','#bc6825','#597fb1','#8a7967','#4c9d82','#553b82','#ba4666','#485e22','#258bba','#c09231']
PREFERRED={'F11xG14':BLUE,'F11xHtL':BLUE,'G14xHtL':BLUE,
 'EtL_polexcalE':ORANGE,'HtLxEtL_pole':ORANGE,
 'EbarTxHtL':GREEN,'HTxcalE':GREEN,'EtL_npxEtL_pole':GREEN,
 'HTxEtL_pole':RED,'HTxHT':RED,
 'F11xF14':PURPLE,'HTxEtL_np':PURPLE}

def configure_styles(theory):
 """Keep a product's style fixed across the four xB bins and both charges.

π+ determines the priority for the five strongest colors. Other products get
separate colors, rather than the former 51-product global dash registry.
"""
 registry={}
 observables=sorted({o for _,o in theory})
 for obs in observables:
  keys=[k for k in theory if k[1]==obs]
  terms={t for key in keys for b in theory[key] for t in b['terms'] if t!='other'}
  reference=theory.get(('pip',obs),theory[keys[0]])
  def size(t):return max((np.max(np.abs(b['arrays'][t])) for b in reference if t in b['arrays']),default=0.)
  order=sorted(terms,key=lambda t:(-size(t),t));styles={};used=set()
  # Reserve preferred colors for the five dominant products in π+.
  for t in order[:5]:
   c=PREFERRED.get(t)
   if c is not None and c not in used:styles[t]=(c,'solid');used.add(c)
  for t in order:
   if t in styles:continue
   color=next((c for c in COLORS if c not in used),None)
   if color is None:raise ValueError('Need more distinct CFF colors for '+obs)
   styles[t]=(color,'solid');used.add(color)
  styles['other']=('#6b6b6b',(0,(4,2)))
  styles['remainder']=('#6b6b6b',(0,(4,2)))
  registry[obs]=styles
  for key in keys:
   for block in theory[key]:block['styles']=styles
 return registry
