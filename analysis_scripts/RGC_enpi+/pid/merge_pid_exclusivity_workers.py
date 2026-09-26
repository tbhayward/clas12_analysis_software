#!/usr/bin/env python3
"""Merge whitespace worker skims from pid_exclusivity_skim.groovy into one ROOT TTree."""
from __future__ import annotations
import argparse, glob, os
import numpy as np
import uproot

COLS = """run event period assigned_pid charge status fd beta chi2pid p px py pz theta phi vz e_p e_px e_py e_pz e_theta e_phi e_vz Ebeam Q2 W xB y t tmin minus_tprime Mx2 analysis_bin mx2_mu mx2_sigma mx2_nsigma pass_qa pass_fd pass_p pass_chi2pid pass_Q2 pass_W pass_y pass_phase_space pass_mx2 pass_final""".split()
PERIOD_CODE = {"Su22":0, "Fa22":1, "Sp23":2}
INTS = {"run","event","assigned_pid","charge","status","fd","analysis_bin","pass_qa","pass_fd","pass_p","pass_chi2pid","pass_Q2","pass_W","pass_y","pass_phase_space","pass_mx2","pass_final"}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("worker_glob")
    ap.add_argument("output_root")
    args=ap.parse_args()
    files=sorted(glob.glob(args.worker_glob))
    if not files: raise SystemExit(f"No workers match {args.worker_glob}")
    os.makedirs(os.path.dirname(os.path.abspath(args.output_root)),exist_ok=True)
    tree=None; total=0
    with uproot.recreate(args.output_root) as fout:
        for fn in files:
            rows=[]
            with open(fn) as f:
                for line in f:
                    if not line.strip() or line.startswith('#'): continue
                    rows.append(line.split())
            if not rows: continue
            a=np.asarray(rows,dtype=object)
            out={}
            for j,c in enumerate(COLS):
                if c=="period": out["period_code"]=np.array([PERIOD_CODE[x] for x in a[:,j]],dtype=np.int8)
                elif c in INTS: out[c]=a[:,j].astype(np.int64)
                else: out[c]=a[:,j].astype(np.float64)
            if tree is None:
                fout["pid_exclusivity"] = out
                tree=fout["pid_exclusivity"]
            else:
                tree.extend(out)
            total += len(a)
            print(f"merged {fn}: {len(a):,} rows")
    print(f"Wrote {total:,} candidates to {args.output_root}:pid_exclusivity")

if __name__ == "__main__": main()
