#!/usr/bin/env python3
"""Launch a t test in a private working directory using one shared build."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--merge-low-t', action='store_true')
    group.add_argument('--photon-t', action='store_true')
    parser.add_argument('--run-dir', required=True, type=Path,
                        help='New, unique directory; existing directories are refused')
    parser.add_argument('--skip-systematics', action='store_true', required=True)
    parser.add_argument('--plot', action='store_true', help='Plot final results after extraction')
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    run = args.run_dir.expanduser().resolve()
    for name in ('dvcs_analysis', 'main_systematics'):
        if not (source / name).is_file():
            parser.error(f'Build first: missing {source / name}')
    # Atomic mkdir prevents two jobs from accidentally sharing a destination.
    run.mkdir(parents=True, exist_ok=False)
    # Snapshot only small code/resources and executable files. ROOT inputs remain
    # at the loader's shared absolute paths; imports are shared read-only inputs.
    for path in source.iterdir():
        if path.is_file() and (path.suffix in ('.py', '.cpp', '.h', '.json', '.sh')
                               or path.name in ('dvcs_analysis', 'main_systematics')):
            shutil.copy2(path, run / path.name)
    for name in ('external_scripts', 'external_pi0_normalization'):
        shutil.copytree(source / name, run / name,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc', 'output'))
    (run / 'imports').symlink_to(source / 'imports', target_is_directory=True)
    mode = '--merge-low-t' if args.merge_low_t else '--photon-t'
    (run / '.dvcs_t_test_run').write_text(json.dumps({
        'source': str(source), 'mode': mode, 'systematics': False,
        'merged_corrections': 'upper valid row, lower-only fallback',
    }, indent=2) + '\n')
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['MPLCONFIGDIR'] = str(run / 'output/.matplotlib')
    # ROOT/C++ stages already allocate up to ten workers. Prevent numerical
    # libraries inside each Python worker from creating another CPU pool.
    for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[name] = '1'
    env['OMP_NUM_THREADS'] = '10'
    print(f'Private run directory: {run}', flush=True)
    result = subprocess.run([str(run / 'dvcs_analysis'), mode, '--skip-systematics'],
                            cwd=run, env=env)
    if result.returncode:
        return result.returncode
    if args.plot:
        return subprocess.run([
            sys.executable, str(run / 'external_scripts/plot_pass2_final_results_by_t.py'),
            str(run / 'output/csvs/dvcs_pass2_analysis.csv'),
            '--output-dir', str(run / 'output/final_results_by_t'),
            '--workers', '10', '--no-pass1', '--reference-kinematics', 'means',
        ], cwd=run, env=env).returncode
    return 0


if __name__ == '__main__':
    sys.exit(main())
