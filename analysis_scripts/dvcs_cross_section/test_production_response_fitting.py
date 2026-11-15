#!/usr/bin/env python3
"""Native original-event cache and training-only fresh optimizer checks."""
import json,tempfile,subprocess,shutil
from pathlib import Path
from unittest.mock import patch
import numpy as np
import fit_production_particle_response as fit
from production_response_stream import stream_fit_cache
import test_production_particle_response as test
with tempfile.TemporaryDirectory() as td:
 d=Path(td);worker,rows,globals_=test.native_build(d)
 for control,name in ((False,'signal'),(True,'control'),(True,'control_swapped')):
  source=d/(name+'.root');test.write_tree(source,rows[name],n=120)
  outputs=[]
  for batch in (17,43):
   job=dict(source=str(source),target=str(d/(name+str(batch)+'.bin')),period='fa18_inb',control=control,beam_energy=10.604,seed=1729,global_cuts=globals_,cache_probability=1.)
   with (d/'log').open('w') as log:m=stream_fit_cache(job,worker,log,batch)
   assert m['entries']==120 and m['records']==120
   raw=Path(job['target']).read_bytes();assert len(raw)==120*(280 if control else 200);outputs.append(raw)
  assert outputs[0]==outputs[1]
  if control:
   dt=np.dtype([('v','<f8',24),('key','<u8'),('photons','<f8',10)]);a=np.frombuffer(outputs[0],dtype=dt);assert np.all(a['photons'][:,9]==1) and np.all(a['v'][:,4]==rows[name]['detector_gamma1'])
 print('PASS: fresh original-event fitting caches, both photon tags/recovery, no event filtering and batch invariance')
 class Engine:
  spec={'periods':['fa18_inb']}
  def evaluate(self,pars,label):
   p=[tuple([0.]*20)]+list(dict.fromkeys(tuple(p) for p in pars if any(p)));scores=[dict(train=sum((v-(.01 if k==3 else 0))**2 for k,v in enumerate(x)),background_multiplier=1.,validation_mc=-sum(v*v for v in x)) for x in p];return p,scores,[],{},None
 grid={k:[0.,.01] if k%4==3 else [0.] for k in range(20)}
 with patch.object(fit.base,'parameter_grid',return_value=grid):best=fit.fit_coordinates(Engine(),d,coarse=1,refine=1)
 assert best[3]==.01 and sum(abs(v) for k,v in enumerate(best) if k!=3)==0
 assert fit.read(d/'candidate_parameters.json')['initial_parameters']=='zero; no imported fit'
 print('PASS: production optimizer starts from zero, preserves shared particle parameters and chooses training score only')
