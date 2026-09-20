#!/usr/bin/env python3
"""Single entry point for the cleaned high-luminosity DVCS workshop projection."""
from pathlib import Path
import argparse, subprocess, sys, shutil, time

def run(cmd):
    print("\n"+"="*100); print("[run] "+" ".join(map(str,cmd))); print("="*100,flush=True)
    subprocess.run(list(map(str,cmd)),check=True)

def main():
    here=Path(__file__).resolve().parent
    p=argparse.ArgumentParser()
    p.add_argument("--replicas",type=int,default=100000)
    p.add_argument("--workers",type=int,default=8)
    p.add_argument("--summary-sample",type=int,default=100000)
    p.add_argument("--scatter-sample",type=int,default=20000)
    p.add_argument("--save-replica-sample",type=int,default=0)
    p.add_argument("--skip-stage1",action="store_true")
    p.add_argument("--skip-cff",action="store_true")
    p.add_argument("--force-derivatives",action="store_true")
    p.add_argument("--clean-output",action="store_true")
    a=p.parse_args(); out=here/"output"; py=sys.executable; t=time.time()
    if a.clean_output and out.exists(): shutil.rmtree(out)
    if not a.skip_stage1:
        run([py,here/"high_luminosity_dvcs_reach.py","--output",out/"stage1"])
    run([py,here/"high_luminosity_dvcs_stage2_pseudodata_pass2.py","--outdir",out/"stage2"])
    if not a.skip_cff:
        run([py,here/"high_luminosity_dvcs_stage2_cff_sensitivity_pass2.py",
             "--input-dir",out/"stage2"/"tables","--outdir",out/"stage2_cff"])
    c=[py,here/"high_luminosity_dvcs_stage3c_volker_dterm.py","--input-dir",out/"stage2"/"tables",
       "--outdir",out/"stage3c","--workers",a.workers]
    if a.force_derivatives: c.append("--force-derivatives")
    run(c)
    run([py,here/"high_luminosity_dvcs_stage3e_nonlinear_replicas.py",
         "--stage2",out/"stage2"/"tables","--stage3c",out/"stage3c"/"tables","--outdir",out/"stage3e",
         "--replicas",a.replicas,"--workers",a.workers,"--summary-sample",a.summary_sample,
         "--scatter-sample",a.scatter_sample,"--save-replica-sample",a.save_replica_sample])
    print(f"\n[done] {out}  elapsed={(time.time()-t)/60:.1f} min")

if __name__=="__main__": main()
