#!/usr/bin/env python3
"""Focused direct-pi0 Mx2 diagnostic; run beside plot_exclusivity_data_dvcs_pi0_mc.py.
Uses the production loader and its preselection, with no modification to cuts.
"""
import argparse
import importlib.util
from pathlib import Path
import sys
import numpy as np


def describe(label, values, xmin=-0.03, xmax=0.03):
    values = np.asarray(values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    inside = finite[(finite >= xmin) & (finite <= xmax)]
    print(f"{label}: selected={len(values):,} finite={len(finite):,} in_Mx2_window={len(inside):,}", flush=True)
    if len(finite):
        quantiles = np.percentile(finite, [0, 1, 5, 25, 50, 75, 95, 99, 100])
        print("  Mx2 quantiles (GeV^2) [0,1,5,25,50,75,95,99,100]%:",
              " ".join(f"{v:.6g}" for v in quantiles), flush=True)
        print(f"  Mx2 near zero |Mx2|<0.03: {len(inside)/len(finite):.5%}", flush=True)
    #endif


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--script', default='plot_exclusivity_data_dvcs_pi0_mc.py')
    parser.add_argument('--global-cuts-json', default='output/jsons/global_cuts_config.json')
    parser.add_argument('--step-size', default='100 MB')
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('dvcs_exclusivity_diag', Path(args.script).resolve())
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    mod.ACTIVE_GLOBAL_CUTS = mod.load_global_cuts_config(args.global_cuts_json)
    variable = next(v for v in mod.PI0_SHAPE_VARIABLES if v.branch == 'Mx2')
    topo = next(t for t in mod.TOPOLOGIES if t.key == 'FD_FD')
    for period in mod.PERIODS:
        if period.key not in ('sp18_inb', 'sp18_out'):
            continue
        #endif
        print(f'\n=== {period.label} direct-pi0 FD-FD ===', flush=True)
        for sample, path in [('data', period.eppi0_data_file), ('MC', period.eppi0_mc_file)]:
            print(f'{sample} input: {path}', flush=True)
            try:
                print('  resolved Mx2:', mod.resolve_variable_branches(path, [variable]), flush=True)
                result = mod.load_selected_arrays_for_file(path, period.key, [topo], args.step_size, [variable])
                describe(sample, result[topo.key]['Mx2'], variable.xmin, variable.xmax)
            except Exception as exc:
                print(f'  ERROR {type(exc).__name__}: {exc}', flush=True)
                raise
            #endtry
        #endfor
    #endfor
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
