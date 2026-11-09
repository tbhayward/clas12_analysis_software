#!/usr/bin/env python3
"""Run data, reconstructed-MC, and paired-MC truth |t| diagnostics in one command."""
import argparse
from pathlib import Path
import subprocess
import sys
import time

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--period",default="fa18_inb",
                    choices=("fa18_inb","fa18_out","sp18_inb","sp18_out","sp19_inb"))
    ap.add_argument("--input-root",type=Path,
                    default=Path("/volatile/clas12/thayward/dvcs_pass2_full_rebuild"))
    ap.add_argument("--output-dir",type=Path,default=Path("output/t_definition_quick_test"))
    ap.add_argument("--max-files",type=int,default=0,
                    help="Files per sample; zero processes all files")
    ap.add_argument("--workers",type=int,default=10)
    ap.add_argument("--cuts-json",type=Path,
                    default=Path(__file__).resolve().parent/"output/jsons/combined_cuts.json")
    args=ap.parse_args()
    if not 1<=args.workers<=10:
        ap.error("--workers must be 1–10")
    scripts=Path(__file__).resolve().parent
    common=["--period",args.period,"--input-root",str(args.input_root),
            "--output-dir",str(args.output_dir),"--max-files",str(args.max_files),
            "--workers",str(args.workers),"--cuts-json",str(args.cuts_json)]
    tasks=(("data",scripts/"compare_t_definitions.py",["--sample","data"]),
           ("reconstructed MC",scripts/"compare_t_definitions.py",["--sample","mc"]),
           ("paired MC truth",scripts/"compare_t_definitions_truth.py",[]))
    for index,(name,script,extra) in enumerate(tasks,1):
        start=time.monotonic()
        print(f"[STAGE {index}/3] {name}: starting",flush=True)
        subprocess.run([sys.executable,str(script),*common,*extra],check=True)
        print(f"[STAGE {index}/3] {name}: done in {(time.monotonic()-start)/60:.1f} min",flush=True)
    print(f"[OUTPUT] {args.output_dir/args.period}/t_definition_counts.png",flush=True)
    print(f"[OUTPUT] {args.output_dir/args.period}/truth_diagnostics.png",flush=True)
    print(f"[OUTPUT] {args.output_dir/args.period}/truth_reconstruction_errors.png",flush=True)

if __name__=="__main__":
    main()
