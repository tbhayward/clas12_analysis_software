#!/usr/bin/env python3
"""Repair rec_detector2 in paired AAOgen MC ROOT files, without reprocessing HIPO.

Creates verified staged ROOT files, retaining original files untouched unless --replace
is explicitly specified. Requires PyROOT on ifarm. Uses ROOT fast basket cloning
for the 159 unaffected branches and creates a corrected detector branch.
"""
import argparse
import glob
import os
import pathlib
import sys
import uuid
from array import array

DEFAULT_ROOT = '/volatile/clas12/thayward/dvcs_pass2_full_rebuild/mc'
PERIODS = ('sp18_inb', 'sp18_out', 'fa18_inb', 'fa18_out', 'sp19_inb')


def inspect(path, ROOT):
    f = ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie():
        raise RuntimeError(f'Cannot open {path}')
    t = f.Get('MCEvents')
    if not t:
        f.Close()
        raise RuntimeError(f'Missing MCEvents in {path}')
    required = ('rec_detector2', 'rec_detector_gamma1', 'rec_kinematics_valid')
    for name in required:
        if not t.GetBranch(name):
            f.Close()
            raise RuntimeError(f'Missing {name} in {path}')
    return f, t


def patch_one(src, dst, ROOT, replace=False, verify_only=False, chunk=1000000):
    f, t = inspect(src, ROOT)
    n = int(t.GetEntries())
    print(f'\nSOURCE {src}\n  entries={n:,}', flush=True)
    if verify_only:
        t.SetBranchStatus('*', 0)
        for b in ('rec_detector2', 'rec_detector_gamma1', 'rec_kinematics_valid'):
            t.SetBranchStatus(b, 1)
        n_valid = n_match = 0
        for i in range(n):
            t.GetEntry(i)
            if t.rec_kinematics_valid >= 0.5:
                n_valid += 1
                n_match += int(t.rec_detector2 == t.rec_detector_gamma1)
            #endif
        #endfor
        print(f'  valid={n_valid:,} detector2==gamma1={n_match:,}', flush=True)
        f.Close()
        return
    #endif
    if dst.exists():
        f.Close()
        raise FileExistsError(f'Refusing to overwrite staged file: {dst}')
    #endif
    dst.parent.mkdir(parents=True, exist_ok=True)
    temp = dst.with_name(dst.name + '.partial-' + uuid.uuid4().hex[:8])
    try:
        out = ROOT.TFile.Open(str(temp), 'RECREATE')
        if not out or out.IsZombie():
            raise RuntimeError(f'Cannot create {temp}')
        #endif
        # Fast clone copies compressed baskets for every unaffected branch.
        t.SetBranchStatus('rec_detector2', 0)
        out.cd()
        new = t.CloneTree(-1, 'fast')
        if not new:
            raise RuntimeError('CloneTree failed')
        #endif
        if int(new.GetEntries()) != n:
            raise RuntimeError('Cloned entry count differs from source')
        #endif
        if new.GetBranch('rec_detector2'):
            raise RuntimeError('Old rec_detector2 unexpectedly survived cloning')
        #endif
        t.SetBranchStatus('*', 0)
        t.SetBranchStatus('rec_detector_gamma1', 1)
        t.SetBranchStatus('rec_kinematics_valid', 1)
        # Detector fields are stored as doubles in the paired ROOT schema.
        gamma = array('d', [0.0]); valid = array('d', [0.0])
        t.SetBranchAddress('rec_detector_gamma1', gamma)
        t.SetBranchAddress('rec_kinematics_valid', valid)
        repaired = array('d', [-1.0])
        branch = new.Branch('rec_detector2', repaired, 'rec_detector2/D')
        n_valid = n_fixed = 0
        for i in range(n):
            t.GetEntry(i)
            if valid[0] >= 0.5:
                repaired[0] = gamma[0]
                n_valid += 1
                n_fixed += int(repaired[0] in (0.0, 1.0))
            else:
                repaired[0] = -1.0
            #endif
            branch.Fill()
            if (i + 1) % chunk == 0:
                print(f'  detector branch filled: {i + 1:,}/{n:,}', flush=True)
            #endif
        #endfor
        if int(branch.GetEntries()) != n:
            raise RuntimeError('New branch entry count differs from tree')
        #endif
        new.Write('', ROOT.TObject.kOverwrite)
        out.Close()
        f.Close()
        # Check the staged file before promoting it.
        check = ROOT.TFile.Open(str(temp), 'READ')
        ct = check.Get('MCEvents') if check else None
        if not ct or int(ct.GetEntries()) != n or not ct.GetBranch('rec_detector2'):
            raise RuntimeError('Staged ROOT verification failed')
        #endif
        if int(ct.GetBranch('rec_detector2').GetEntries()) != n:
            raise RuntimeError('Staged detector branch has incorrect entry count')
        #endif
        for i in sorted(set((0, n // 2, n - 1))):
            if i >= 0 and n:
                ct.GetEntry(i)
                expected = float(ct.rec_detector_gamma1) if ct.rec_kinematics_valid >= 0.5 else -1.0
                if float(ct.rec_detector2) != expected:
                    raise RuntimeError(f'Verification mismatch at entry {i}')
                #endif
            #endif
        #endfor
        check.Close()
        os.replace(temp, dst)
        print(f'  VERIFIED staged={dst}\n  valid={n_valid:,} valid photon IDs={n_fixed:,}', flush=True)
        if replace:
            backup = src.with_name(src.name + '.pre_detector2_fix')
            if backup.exists():
                raise FileExistsError(f'Backup already exists: {backup}')
            #endif
            os.replace(src, backup)
            try:
                os.replace(dst, src)
            except Exception:
                os.replace(backup, src)
                raise
            #endtry
            print(f'  INSTALLED {src}; original backup={backup}', flush=True)
        #endif
    finally:
        f.Close()
        if temp.exists():
            print(f'  Removing incomplete staging file {temp}', flush=True)
            temp.unlink()
        #endif
    #endtry


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', default=DEFAULT_ROOT)
    p.add_argument('--periods', nargs='+', choices=PERIODS, default=list(PERIODS))
    p.add_argument('--max-files', type=int, default=1, help='Per period; 0 means all')
    p.add_argument('--stage-dir', default=None, help='Destination root for corrected files (default: sibling staging directory)')
    p.add_argument('--replace', action='store_true', help='Install verified copies and retain .pre_detector2_fix originals (requires disk space)')
    p.add_argument('--verify-only', action='store_true', help='Check existing detector2 values without modifying files')
    args = p.parse_args()
    try:
        import ROOT
    except ImportError as exc:
        raise SystemExit(f'PyROOT unavailable; run on ifarm with ROOT environment: {exc}')
    #endtry
    ROOT.gROOT.SetBatch(True)
    root = pathlib.Path(args.root)
    stage = pathlib.Path(args.stage_dir) if args.stage_dir else root.parent / 'pi0_mc_detector2_staging'
    for period in args.periods:
        pattern = str(root / period / 'pi0_mc' / 'nominal' / 'parts' / '*.root')
        paths = sorted(pathlib.Path(s) for s in glob.glob(pattern) if 'nobkg' not in s.lower())
        if args.max_files:
            paths = paths[:args.max_files]
        #endif
        print(f'PERIOD {period}: selected {len(paths)} files', flush=True)
        for src in paths:
            dst = stage / period / src.name
            patch_one(src, dst, ROOT, replace=args.replace, verify_only=args.verify_only)
        #endfor
    #endfor

if __name__ == '__main__':
    main()
#endif
