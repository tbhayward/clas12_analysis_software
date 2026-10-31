#!/usr/bin/env python3
"""One invocation per ifarm node: run all configured periods for validated paired MC and data modes.
After both nodes finish, node 1 optionally validates/merges both outputs. Never modifies QADB or cuts.
Includes pi0_as_dvcs on the AAOgen input, subject to pilot validation.
"""
import argparse, collections, json, pathlib, subprocess, sys, time, shutil
from campaign_manifest import inventory, freeze, verify
P=pathlib.Path
PERIODS=('sp18_inb','sp18_out','fa18_inb','fa18_out','sp19_inb')
def call(cmd):
 print('[COMMAND]', ' '.join(map(str,cmd)),flush=True)
 subprocess.run(list(map(str,cmd)),check=True)
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--node',type=int,required=True,choices=(1,2));p.add_argument('--workers',type=int,default=10)
 p.add_argument('--output-root',default='rebuild_full_20261028');p.add_argument('--max-events',type=int,default=0)
 p.add_argument('--execute',action='store_true');p.add_argument('--merge-on-node1',action='store_true',help='Node 1 waits for both nodes then merges all periods and MC groups')
 p.add_argument('--wait-hours',type=float,default=72)
 a=p.parse_args()
 if a.workers<1 or a.max_events<0:p.error('workers must be positive and max-events nonnegative')
 root=P.cwd();base=root/a.output_root;cfg=json.loads((root/'processing_scripts/campaign_config.json').read_text())
 manifest=inventory(cfg,a.max_events)
 jobs=[j for j in cfg['jobs'] if j['period'] in PERIODS and j['mode'] in ('dvcs_data','pi0_data','dvcs_rec','pi0_rec')]
 missing=[str(j['input']) for j in jobs if not P(j['input']).is_dir() or not any(P(j['input']).rglob('*.hipo'))]
 print('[PERIOD_ENERGIES]',sorted(set((j['period'],j['energy']) for j in jobs)),flush=True)
 print('[MODES] dvcs_data pi0_data dvcs_mc pi0_mc pi0_as_dvcs (rec-only)' ,flush=True)
 print('[INPUTS]',len(jobs),'configured directories; missing/empty',len(missing),flush=True)
 for x in missing:print('[MISSING]',x,flush=True)
 if missing:raise SystemExit('Refusing incomplete campaign; fix input paths first')
 if not (root/'processing_scripts/rebuild_pilot/campaign_worker.py').is_file():p.error('Missing MC worker')
 required_sources=(
  'data_pilot_worker.py',
  'pi0_as_dvcs_worker.py',
  'finalize_data_campaign.py',
  'convert_data_stream.cpp',
  'validate_reconstructed_root.cpp',
  'reconstructed_data_schema.h')
 for name in required_sources:
  if not (root/'processing_scripts/rebuild_pilot'/name).is_file():p.error('Missing '+name)
 if not a.execute:
  print('[DRY RUN] no processing performed; add --execute on each node',flush=True);return
 if base.exists() and not base.is_dir():p.error('output root is not a directory')
 base.mkdir(parents=True,exist_ok=True)
 manifest_hash=freeze(base,manifest)
 # Each worker has its own output root, preventing collisions and allowing safe independent restarts.
 mc=base/'mc';data=base/'data';leakage=base/'leakage'
 commands=[
  [sys.executable,'processing_scripts/rebuild_pilot/campaign_worker.py','--node',a.node,'--nodes',2,'--workers',a.workers,'--max-events',a.max_events,'--output-root',mc,'--execute'],
  [sys.executable,'processing_scripts/rebuild_pilot/data_pilot_worker.py','--node',a.node,'--nodes',2,'--workers',a.workers,'--max-events',a.max_events,'--output-root',data,'--periods',*PERIODS,'--execute','--retry-stale'],
  [sys.executable,'processing_scripts/rebuild_pilot/pi0_as_dvcs_worker.py','--node',a.node,'--nodes',2,'--workers',a.workers,'--max-events',a.max_events,'--output-root',leakage,'--execute','--retry-stale']]
 for command in commands:
  if inventory(cfg,a.max_events)!=manifest:raise SystemExit('Inputs changed since manifest freeze')
  command.extend(['--periods',*PERIODS]) if '--periods' not in command else None
  call(command)
 verify(base,manifest,node=a.node)
 done=base/f'.node{a.node}.complete.json'
 done.write_text(json.dumps({'node':a.node,'manifest_sha256':manifest_hash,'max_events':a.max_events,'periods':PERIODS,'modes':['dvcs_mc','pi0_mc','pi0_as_dvcs','dvcs_data','pi0_data']})+'\n')
 print('[NODE_COMPLETE]',a.node,flush=True)
 if a.node!=1 or not a.merge_on_node1:return
 deadline=time.monotonic()+a.wait_hours*3600
 while True:
  markers=[]
  for n in (1,2):
   f=base/f'.node{n}.complete.json'
   if f.is_file():
    try:
     m=json.loads(f.read_text())
     if m.get('node')==n and m.get('max_events')==a.max_events and m.get('manifest_sha256')==manifest_hash:markers.append(n)
    except (ValueError,OSError):pass
  if len(markers)==2:break
  if time.monotonic()>deadline:raise SystemExit('Timed out waiting for second node; processing outputs preserved')
  print('[WAIT] node 1 waiting for node 2 to complete',flush=True);time.sleep(120)
 verify(base,manifest)
 if inventory(cfg,a.max_events)!=manifest:raise SystemExit('Inputs changed since manifest freeze')
 # Counts come from the frozen inventory.
 for period in PERIODS:
  expected=sum(r['period']==period and r['mode'] in ('dvcs_data','pi0_data') for r in manifest['records'])
  call([sys.executable,'processing_scripts/rebuild_pilot/finalize_data_campaign.py','--output-root',data,'--period',period,'--expected-files',expected,'--max-events',a.max_events,'--merge'])
 for period in PERIODS:
  expected=sum(r['period']==period and r['mode']=='pi0_as_dvcs' for r in manifest['records'])
  call([sys.executable,'processing_scripts/rebuild_pilot/finalize_data_campaign.py','--output-root',leakage,'--period',period,'--expected-files',expected,'--max-events',a.max_events,'--only','pi0_as_dvcs','--merge'])
 call([sys.executable,'processing_scripts/rebuild_pilot/finalize_mc_campaign.py','--output-root',mc,'--merge']+(['--allow-pilot'] if a.max_events else []))
 print('[ALL_COMPLETE] all configured periods processed and merged for five modes',flush=True)
if __name__=='__main__':main()
