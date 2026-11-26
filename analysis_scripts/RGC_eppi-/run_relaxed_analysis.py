#!/usr/bin/env python3
"""Rerun the existing 24-bin polarized analysis on the relaxed campaign."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGES = ('selection', 'audit', 'kinematics', 'dilution', 'asymmetry')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--datasets', type=Path, default=HERE/'imports/datasets.json')
    parser.add_argument('--input-root', type=Path)
    parser.add_argument('--output-root', type=Path, default=HERE/'output_loose_20261124')
    parser.add_argument('--run-info', type=Path, default=HERE/'imports/clas12_run_info_ND3_source.csv')
    parser.add_argument('--review-output', type=Path, default=HERE/'review_loose_20261124', help='Compact review copy and sibling tar.gz, updated after each successful stage')
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--stage', choices=('all',)+STAGES, default='all')
    parser.add_argument('--status', action='store_true', help='Show readable selection progress without running analysis')
    parser.add_argument('--dry-run', action='store_true', help='Print exact commands without reading inputs or running stages')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be positive')
    config = json.loads(args.datasets.read_text())
    source = args.input_root or Path(config['campaign'])/'root'
    out = args.output_root.resolve()
    selected = out/'channel_selection/selected_inputs.json'
    dilution = out/'dilution_factor/joint_period_method1.json'
    selection = [sys.executable, str(HERE/'channel_selection/run_channel_selection.py'), '--datasets', str(args.datasets.resolve()), '--input-root', str(source.resolve()), '--output', str(out/'channel_selection')]
    if args.status:
        subprocess.run(selection+['--status'], check=True)
        return
    commands = {
        'selection': selection+['--workers', str(args.workers)],
        'audit': selection+['--audit'],
        'kinematics': [sys.executable, str(HERE/'kinematic_distributions/plot_kinematic_distributions.py'), '--inputs', str(selected), '--output', str(out/'kinematic_distributions'), '--workers', str(args.workers)],
        'dilution': [sys.executable, str(HERE/'dilution_factor/determine_dilution_factor.py'), '--selected-inputs', str(selected), '--run-info', str(args.run_info.resolve()), '--output', str(out/'dilution_factor'), '--joint-workers', str(args.workers)],
        'asymmetry': [sys.executable, str(HERE/'asymmetry_extraction/extract_structure_function_ratios.py'), '--selected-inputs', str(selected), '--run-info', str(args.run_info.resolve()), '--dilution', str(dilution), '--output', str(out/'asymmetry_extraction'), '--workers', str(args.workers)],
    }
    stages = STAGES if args.stage == 'all' else (args.stage,)
    if not args.dry_run and 'selection' in stages:
        missing = [str(source/name) for name in config['files'].values() if not (source/name).is_file()]
        if missing:
            parser.error('Wait for all 19 final merges. Missing inputs:\n'+'\n'.join(missing))
    env = os.environ.copy()
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[name] = '1'
    env['MPLBACKEND'] = 'Agg'
    import shlex
    for stage in stages:
        print('\n'+stage.upper()+': '+shlex.join(commands[stage]), flush=True)
        if not args.dry_run:
            subprocess.run(commands[stage], check=True, env=env)
            from export_review import export_review
            export_review(out, args.review_output, args.run_info, args.datasets)
    if not args.dry_run:
        print('\nResults: '+str(out), flush=True)


if __name__ == '__main__':
    main()
