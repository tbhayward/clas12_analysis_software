#!/usr/bin/env python3
"""Repair paired AAOgen MC rec_detector2 using a full fast clone.

Unlike v1, NEVER disable a source branch before CloneTree('fast'): that
combination can produce an invalid TTree in ROOT 6.36/PyROOT. Clone all
branches, rename the old detector branch, then append the corrected branch.
Original files are untouched unless --replace is explicitly requested.
"""
import argparse
import glob
import os
from pathlib import Path
from array import array
from uuid import uuid4

PERIODS = ('sp18_inb', 'sp18_out', 'fa18_inb', 'fa18_out', 'sp19_inb')
ROOT_DIR = Path('/volatile/clas12/thayward/dvcs_pass2_full_rebuild/mc')


def patch(src, dst, ROOT, replace=False, progress=5000000):
    print(f'\nSOURCE {src}', flush=True)
    fin = ROOT.TFile.Open(str(src), 'READ')
    if not fin or fin.IsZombie():
        raise RuntimeError('Cannot open source ROOT file')
    #endif
    tree = fin.Get('MCEvents')
    if not tree:
        raise RuntimeError('Missing MCEvents')
    #endif
    for name in ('rec_detector2', 'rec_detector_gamma1', 'rec_kinematics_valid'):
        if not tree.GetBranch(name):
            raise RuntimeError(f'Missing branch {name}')
        #endif
    #endfor
    n = int(tree.GetEntries())
    print(f'  entries={n:,}', flush=True)
    if dst.exists():
        raise FileExistsError(f'Staging destination already exists: {dst}')
    #endif
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + '.partial-' + uuid4().hex[:8])
    fout = None
    try:
        fout = ROOT.TFile.Open(str(tmp), 'RECREATE')
        if not fout or fout.IsZombie():
            raise RuntimeError('Cannot create output ROOT file')
        #endif
        # Crucial difference from v1: clone ALL source branches with no
        # disabled branch. ROOT fast cloning copies compressed baskets.
        tree.SetBranchStatus('*', 1)
        fout.cd()
        cloned = tree.CloneTree(-1, 'fast')
        if not cloned:
            raise RuntimeError('Fast CloneTree failed')
        #endif
        if int(cloned.GetEntries()) != n:
            raise RuntimeError('Fast clone entry count mismatch')
        #endif
        print('  fast clone completed', flush=True)
        # Retain the original values for audit, under a new name.
        old = cloned.GetBranch('rec_detector2')
        if not old:
            raise RuntimeError('Original detector2 branch missing after full clone')
        #endif
        old.SetName('rec_detector2_original')
        # Rename leaf too, so there is no ambiguous leaf name.
        leaf = old.GetLeaf('rec_detector2')
        if leaf:
            leaf.SetName('rec_detector2_original')
        #endif
        if cloned.GetBranch('rec_detector2'):
            raise RuntimeError('Branch rename did not remove old name')
        #endif
        corrected = array('d', [-1.0])
        new_branch = cloned.Branch('rec_detector2', corrected, 'rec_detector2/D')
        if not new_branch:
            raise RuntimeError('Cannot create replacement detector branch')
        #endif
        # Only read two small source branches for the repair loop.
        tree.SetBranchStatus('*', 0)
        tree.SetBranchStatus('rec_detector_gamma1', 1)
        tree.SetBranchStatus('rec_kinematics_valid', 1)
        gamma = array('d', [-1.0])
        valid = array('d', [0.0])
        if tree.SetBranchAddress('rec_detector_gamma1', gamma) < 0:
            raise RuntimeError('Cannot bind rec_detector_gamma1')
        #endif
        if tree.SetBranchAddress('rec_kinematics_valid', valid) < 0:
            raise RuntimeError('Cannot bind rec_kinematics_valid')
        #endif
        nvalid = 0
        nbad = 0
        for i in range(n):
            tree.GetEntry(i)
            if valid[0] >= 0.5:
                corrected[0] = gamma[0]
                nvalid += 1
                if corrected[0] not in (0.0, 1.0):
                    nbad += 1
                #endif
            else:
                corrected[0] = -1.0
            #endif
            new_branch.Fill()
            if progress and (i + 1) % progress == 0:
                print(f'  filled {i+1:,}/{n:,}', flush=True)
            #endif
        #endfor
        if nbad:
            raise RuntimeError(f'{nbad} valid events had invalid gamma1 detector IDs')
        #endif
        if int(new_branch.GetEntries()) != n:
            raise RuntimeError('Replacement branch entry count mismatch')
        #endif
        fout.cd()
        cloned.Write('MCEvents', ROOT.TObject.kOverwrite)
        fout.Close()
        fout = None
        fin.Close()
        print(f'  valid={nvalid:,}; validating staged file', flush=True)
        check = ROOT.TFile.Open(str(tmp), 'READ')
        if not check or check.IsZombie():
            raise RuntimeError('Cannot reopen staged file')
        #endif
        ct = check.Get('MCEvents')
        if not ct or int(ct.GetEntries()) != n:
            raise RuntimeError('Staged MCEvents entry count mismatch')
        #endif
        for name in ('rec_detector2', 'rec_detector2_original', 'rec_detector_gamma1', 'rec_kinematics_valid'):
            if not ct.GetBranch(name):
                raise RuntimeError(f'Staged file missing {name}')
            #endif
        #endfor
        if int(ct.GetBranch('rec_detector2').GetEntries()) != n:
            raise RuntimeError('New branch has wrong entry count')
        #endif
        for i in sorted(set([0, 1, n // 3, n // 2, n - 1])):
            if 0 <= i < n:
                ct.GetEntry(i)
                expected = float(ct.rec_detector_gamma1) if float(ct.rec_kinematics_valid) >= 0.5 else -1.0
                if float(ct.rec_detector2) != expected:
                    raise RuntimeError(f'Value mismatch at entry {i}: {ct.rec_detector2} != {expected}')
                #endif
            #endif
        #endfor
        check.Close()
        os.replace(tmp, dst)
        print(f'  VERIFIED: {dst}', flush=True)
        if replace:
            backup = src.with_name(src.name + '.pre_detector2_fix')
            if backup.exists():
                raise FileExistsError(f'Backup exists: {backup}')
            #endif
            os.replace(src, backup)
            try:
                os.replace(dst, src)
            except Exception:
                os.replace(backup, src)
                raise
            #endtry
            print(f'  INSTALLED: {src}; backup: {backup}', flush=True)
        #endif
    finally:
        if fout:
            fout.Close()
        #endif
        if tmp.exists():
            print(f'  Removing incomplete output {tmp}', flush=True)
            tmp.unlink()
        #endif
    #endtry


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT_DIR)
    parser.add_argument('--periods', nargs='+', choices=PERIODS, default=list(PERIODS))
    parser.add_argument('--max-files', type=int, default=1, help='Per period, 0 means all')
    parser.add_argument('--stage-dir', type=Path)
    parser.add_argument('--replace', action='store_true')
    parser.add_argument('--progress', type=int, default=5000000)
    args = parser.parse_args()
    import ROOT
    ROOT.gROOT.SetBatch(True)
    stage = args.stage_dir or args.root.parent / 'pi0_mc_detector2_staging'
    for period in args.periods:
        pattern = str(args.root / period / 'pi0_mc' / 'nominal' / 'parts' / '*.root')
        paths = [Path(p) for p in sorted(glob.glob(pattern)) if 'nobkg' not in p.lower() and not p.endswith('.pre_detector2_fix')]
        if args.max_files > 0:
            paths = paths[:args.max_files]
        #endif
        print(f'PERIOD {period}: selected {len(paths)} files', flush=True)
        for src in paths:
            patch(src, stage / period / src.name, ROOT, replace=args.replace, progress=args.progress)
        #endfor
    #endfor


if __name__ == '__main__':
    main()
#endif
