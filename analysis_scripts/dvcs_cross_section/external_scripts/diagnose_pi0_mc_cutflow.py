#!/usr/bin/env python3
"""Trace zero reconstructed AAOgen pi0 MC population, using production cuts.
Run from dvcs_cross_section; needs uproot, numpy and production fitter script.
"""
import argparse
import glob
import importlib.util
import os
import warnings
from pathlib import Path
import numpy as np


def load_production(path):
    spec = importlib.util.spec_from_file_location('dvcs_exclusivity_production', path)
    mod = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--periods', nargs='+', default=['sp18_inb', 'sp18_out', 'fa18_inb'])
    ap.add_argument('--max-files', type=int, default=2, help='Files per period (0 = all)')
    ap.add_argument('--max-entries', type=int, default=500000, help='Entries per period (0 = all)')
    ap.add_argument('--step-size', default='50 MB')
    ap.add_argument('--production-script', default='plot_exclusivity_data_dvcs_pi0_mc.py')
    ap.add_argument('--global-cuts-json', default='output/jsons/global_cuts_config.json')
    ap.add_argument('--root', default='/volatile/clas12/thayward/dvcs_pass2_full_rebuild')
    args = ap.parse_args()
    import uproot
    prod = load_production(os.path.abspath(args.production_script))
    cfg = prod.load_global_cuts_config(args.global_cuts_json)
    prod.ACTIVE_GLOBAL_CUTS = cfg
    print('Production global cuts:', args.global_cuts_json, 'tag=', cfg.analysis_tag, flush=True)
    for period in args.periods:
        current = prod.REBUILT_NOMINAL_CURRENT[period]
        pattern = f'{args.root}/mc/{period}/pi0_mc/nominal/parts/*_{current}_*.root'
        all_files = sorted(glob.glob(pattern))
        files = all_files if args.max_files == 0 else all_files[:args.max_files]
        print(f'\n=== {period} nominal pi0 MC ({current}) ===', flush=True)
        print(f'pattern: {pattern}\nmatched={len(all_files)}, inspecting={len(files)}', flush=True)
        if not files:
            print('FATAL: no files matched; investigate filename/current selection', flush=True)
            continue
        for path in files:
            print(' file:', os.path.basename(path), flush=True)
        with uproot.open(files[0]) as rootfile:
            if 'MCEvents' not in rootfile:
                print('FATAL: MCEvents absent; keys=', rootfile.keys(), flush=True)
                continue
            tree = rootfile['MCEvents']
            names = set(tree.keys())
            print(' first file entries=', tree.num_entries, 'branch count=', len(names), flush=True)
            required = sorted(set(prod.global_cut_required_branches() + ['Mx2', 'Q2', 'x', 'W', 'p1_p', 'p1_theta', 'p2_p', 'p2_theta']))
            actual = {name: 'rec_' + name for name in required if 'rec_' + name in names}
            missing = [name for name in required if name not in actual]
            if missing:
                print(' MISSING reconstructed branches:', missing, flush=True)
            print(' validity branch:', 'rec_kinematics_valid' in names, 'gen validity:', 'gen_kinematics_valid' in names, flush=True)
            print(' detector branch names:', sorted(n for n in names if 'detector' in n)[:20], flush=True)
            if 'rec_kinematics_valid' not in names:
                continue
        needed = sorted(set(actual.values()) | {'rec_kinematics_valid'})
        stages = ['all', 'rec_valid', 'FD-FD', 'finite_t1_open', 't1_limit', 'open_angle', 'pTmiss', 'sector_e', 'sector_p', 'sector_gamma', 'distinct_sectors', 'production_global']
        totals = {key: 0 for key in stages}
        independent = {key: 0 for key in stages[2:]}
        seen = 0
        missing_needed = bool(missing)
        if missing_needed:
            print(' Skipping cutflow: required rec_* branches missing.', flush=True)
            continue
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message=r"^'(?:t|t1|t2)' is both an alias and a branch name$")
            for chunk in uproot.iterate([f'{p}:MCEvents' for p in files], expressions=needed, step_size=args.step_size, library='np'):
                n = len(chunk['rec_kinematics_valid'])
                if args.max_entries and seen + n > args.max_entries:
                    n = args.max_entries - seen
                if n <= 0:
                    break
                arr = {name: np.asarray(chunk[phys][:n]) for name, phys in actual.items()}
                valid_raw = np.asarray(chunk['rec_kinematics_valid'][:n])
                valid = np.isfinite(valid_raw) & (valid_raw >= 0.5)
                d1 = arr['detector1']; d2 = arr['detector2']
                topo = (d1 == 1) & (d2 == 1)
                t1 = np.asarray(arr['t1'], dtype=float)
                opening = np.asarray(arr['open_angle_ep2'], dtype=float)
                finite = np.isfinite(t1) & np.isfinite(opening)
                tcut = (-t1) < cfg.t1_abs_max
                ocut = opening > cfg.open_angle_min_deg
                pt = np.isfinite(arr['pTmiss']) & (arr['pTmiss'] <= cfg.pTmiss_max) if cfg.enable_pTmiss_cut else np.ones(n, dtype=bool)
                es = prod.fd_sector_from_phi(arr['e_phi'])
                ps = prod.fd_sector_from_phi(arr['p1_phi'])
                gs = prod.fd_sector_from_phi(arr['p2_phi'])
                eok = es >= 1
                pok = ps >= 1
                gok = gs >= 1
                distinct = (es != ps) & (es != gs) & (ps != gs)
                cuts = [('all', np.ones(n, dtype=bool)), ('rec_valid', valid), ('FD-FD', topo), ('finite_t1_open', finite), ('t1_limit', tcut), ('open_angle', ocut), ('pTmiss', pt), ('sector_e', eok), ('sector_p', pok), ('sector_gamma', gok), ('distinct_sectors', distinct)]
                cumulative = np.ones(n, dtype=bool)
                for label, condition in cuts:
                    cumulative &= condition
                    totals[label] += int(np.count_nonzero(cumulative))
                    if label in independent:
                        independent[label] += int(np.count_nonzero(valid & topo & condition))
                    #endif
                production = prod.global_cut_mask_uproot(arr, period, d1, d2)
                totals['production_global'] += int(np.count_nonzero(valid & topo & production))
                independent['production_global'] += int(np.count_nonzero(valid & topo & production))
                seen += n
                if seen == n:
                    print(' initial validity values:', np.unique(valid_raw, return_counts=True), flush=True)
                    print(' initial detector1 values (valid):', np.unique(d1[valid], return_counts=True), flush=True)
                    print(' initial detector2 values (valid):', np.unique(d2[valid], return_counts=True), flush=True)
                    print(' initial detector pairs (valid):', np.unique(np.stack([d1[valid], d2[valid]], axis=1), axis=0, return_counts=True) if np.any(valid) else 'none', flush=True)
                    if np.any(valid & topo):
                        ix = valid & topo
                        for name in ['t1', 'open_angle_ep2', 'Mx2', 'p1_theta', 'p2_theta']:
                            values = np.asarray(arr[name][ix], dtype=float)
                            finite_values = values[np.isfinite(values)]
                            print(' ', name, 'valid FD-FD quantiles:', np.percentile(finite_values, [0, 1, 50, 99, 100]) if len(finite_values) else 'none', flush=True)
                        #endfor
                    #endif
                #endif
                if args.max_entries and seen >= args.max_entries:
                    break
                #endif
            #endfor
        #endwith
        print(f' Cutflow over {seen:,} paired entries:', flush=True)
        for label in stages:
            print(f'  {label:22s} cumulative={totals[label]:>10,} independent(valid FD-FD)={independent.get(label, 0):>10,}', flush=True)
        #endfor
        if totals['rec_valid'] == 0:
            print(' CONCLUSION: zero valid reconstructed candidates; inspect MC reconstruction and flag production.', flush=True)
        elif totals['FD-FD'] == 0:
            print(' CONCLUSION: valid rec events exist, but none have FD proton + FD photon.', flush=True)
        elif totals['production_global'] == 0:
            print(' CONCLUSION: FD-FD exists but production global cuts reject all; inspect independent cut counts.', flush=True)
        else:
            print(' CONCLUSION: baseline MC population exists; investigate later iterative cuts or previous diagnostic assumptions.', flush=True)
        #endif
    #endfor


if __name__ == '__main__':
    main()
