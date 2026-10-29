#!/usr/bin/env python3
"""Resume all five periods of DVCS and pi0 data for one node; preserves .root.ok outputs."""
import argparse
import pathlib
import subprocess
import sys

PERIODS = ('sp18_inb', 'sp18_out', 'fa18_inb', 'fa18_out', 'sp19_inb')

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--node', type=int, choices=(1, 2), required=True)
    ap.add_argument('--workers', type=int, default=10)
    ap.add_argument('--max-events', type=int, default=0)
    ap.add_argument('--output-root', default='/volatile/clas12/thayward/dvcs_pass2_full_rebuild')
    ap.add_argument('--execute', action='store_true')
    args = ap.parse_args()
    if args.workers < 1 or args.max_events < 0:
        ap.error('workers must be positive and max-events nonnegative')
    worker = pathlib.Path('processing_scripts/rebuild_pilot/data_pilot_worker.py')
    if not worker.is_file():
        ap.error('Run from clas12_analysis_software repository root')
    cmd = [sys.executable, str(worker), '--node', str(args.node), '--nodes', '2',
           '--workers', str(args.workers), '--max-events', str(args.max_events),
           '--output-root', str(pathlib.Path(args.output_root) / 'data'),
           '--periods', *PERIODS]
    if args.execute:
        cmd += ['--execute', '--retry-stale']
    print('[COMMAND]', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True)

if __name__ == '__main__':
    main()
