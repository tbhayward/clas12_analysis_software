#!/usr/bin/env python3
"""Repair existing ROOT branch metadata in-place; never re-read HIPO. Requires root-config and g++."""
import argparse
import subprocess
import tempfile
from pathlib import Path
from schema import SCHEMAS

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',default='/volatile/clas12/thayward/dvcs_pass2_full_rebuild')
    ap.add_argument('--execute',action='store_true',help='Without this flag only print inventory')
    ap.add_argument('--include-leakage',action='store_true',help='Also repair pi0_as_dvcs trees')
    a=ap.parse_args()
    base=Path(a.root)
    modes=[('data','dvcs_data'),('data','pi0_data')]
    if a.include_leakage:modes.append(('leakage','pi0_as_dvcs'))
    jobs=[]
    for section,mode in modes:
        for p in sorted((base/section).glob(f'*/{mode}/parts/*.root')):
            if p.name.endswith('.partial.root') or not Path(str(p)+'.ok').exists():continue
            jobs.append((p,mode))
        #endfor
    #endfor
    print(f'[PLAN] files={len(jobs)} execute={a.execute}',flush=True)
    if not a.execute:return
    here=Path(__file__).resolve().parent
    exe=here/'repair_data_branches'
    flags=subprocess.check_output(['root-config','--cflags','--libs'],text=True).split()
    subprocess.run(['g++','-O2','-std=c++17',str(here/'repair_data_branches.cpp'),'-o',str(exe),*flags],check=True)
    with tempfile.TemporaryDirectory() as td:
        paths={}
        for mode,names in SCHEMAS.items():
            s=Path(td)/(mode+'.txt');s.write_text('\n'.join(names)+'\n');paths[mode]=s
        #endfor
        failures=[]
        for i,(p,mode) in enumerate(jobs,1):
            result=subprocess.run([str(exe),str(p),str(paths[mode])])
            if result.returncode:failures.append(str(p))
            if i%50==0:print(f'[PROGRESS] {i}/{len(jobs)} failed={len(failures)}',flush=True)
        #endfor
    print(f'[FINAL] files={len(jobs)} failed={len(failures)}',flush=True)
    if failures:raise SystemExit('Failed files:\n'+'\n'.join(failures))

if __name__=='__main__':main()
