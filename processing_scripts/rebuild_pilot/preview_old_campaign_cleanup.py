#!/usr/bin/env python3
"""Preview old campaign output trees; optional deletion ONLY of explicitly named paths.
Never deletes sources, configuration, code or a new full-statistics output tree.
"""
import argparse,pathlib,shutil,sys
P=pathlib.Path
DEFAULT=('rebuild_data_pilot_1000','data_root_audit_1000','rebuild_scheduler_test','rebuild_data_pilot_v1','rebuild_data_pilot_1000_backup')
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--paths',nargs='+',default=list(DEFAULT),help='Relative old output directories; review before deletion')
 p.add_argument('--delete',action='store_true',help='Actually delete the explicitly named old output directories')
 a=p.parse_args();root=P.cwd().resolve()
 forbidden={'processing_scripts','processing_classes','coatjava','.git','output','outputs','data','mc','rebuild_full_20261028'}
 for raw in a.paths:
  x=P(raw)
  if x.is_absolute() or len(x.parts)!=1 or x.name in forbidden or x.name.startswith('.'):
   p.error('Unsafe path: '+raw)
  target=root/x
  if not target.exists():print('[ABSENT]',target);continue
  if not target.is_dir() or target.is_symlink():p.error('Not a regular directory: '+raw)
  count=sum(1 for _ in target.rglob('*'))
  print('[DELETE]' if a.delete else '[PREVIEW]',target,'items',count,flush=True)
  if a.delete:shutil.rmtree(target)
 if not a.delete:print('[DRY RUN] Nothing deleted. Review paths, then add --delete.',flush=True)
if __name__=='__main__':main()
