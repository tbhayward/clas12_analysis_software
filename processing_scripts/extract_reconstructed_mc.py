#!/usr/bin/env python3
"""Create sibling reconstructed-only ROOT files from existing paired MC files.

Requires PyROOT. Default is a dry run. Originals are opened read-only.

Standard companion convention (use for existing and future campaigns):
    sample.root -> sample.rec_only.root, tree PhysicsEvents
Run this exporter after a future campaign as well to produce the same interface.
The campaign converter's embedded PhysicsEvents tree is also compatible, but
does not by itself create a separate companion file.
"""
import argparse
from array import array
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
from pathlib import Path
import uuid

SUFFIX = '.rec_only.root'
VERSION = 1


def fingerprint(path):
    stat = path.stat()
    return dict(version=VERSION, source=str(path.resolve()),
                source_size=stat.st_size, source_mtime_ns=stat.st_mtime_ns)


def inspect_output(ROOT, path, identity):
    f = ROOT.TFile.Open(str(path), 'READ')
    if not f or f.IsZombie():
        raise RuntimeError('Cannot open output: ' + str(path))
    try:
        meta = f.Get('reconstructed_export_metadata')
        tree = f.Get('PhysicsEvents')
        if not meta or not tree:
            raise RuntimeError('Missing output tree or provenance: ' + str(path))
        data = json.loads(str(meta.GetTitle()))
        if any(data.get(k) != v for k, v in identity.items()):
            raise RuntimeError('Output belongs to another source/version: ' + str(path))
        if tree.GetEntries() != data['selected_rows']:
            raise RuntimeError('Output entry count mismatch: ' + str(path))
        if [b.GetName() for b in tree.GetListOfBranches()] != data['branches']:
            raise RuntimeError('Output schema mismatch: ' + str(path))
        return data
    finally:
        f.Close()


def extract(source_string):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    source = Path(source_string)
    dest = source.with_name(source.stem + SUFFIX)
    identity = fingerprint(source)
    if dest.exists():
        data = inspect_output(ROOT, dest, identity)
        return f'[SKIP_VALID] {dest} rows={data["selected_rows"]}'
    lock = Path(str(dest) + '.lock')
    # Exclusive creation prevents concurrent exporters from publishing the same file.
    descriptor = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    os.close(descriptor)
    temporary = dest.with_name(dest.name + '.' + uuid.uuid4().hex + '.partial')
    input_file = output_file = None
    try:
        input_file = ROOT.TFile.Open(str(source), 'READ')
        if not input_file or input_file.IsZombie() or input_file.TestBit(ROOT.TFile.kRecovered):
            raise RuntimeError('Unreadable/recovered input: ' + str(source))
        paired = input_file.Get('MCEvents')
        if not paired or not paired.InheritsFrom('TTree'):
            raise RuntimeError('Missing MCEvents: ' + str(source))
        rec_names = [b.GetName() for b in paired.GetListOfBranches()
                     if b.GetName().startswith('rec_') and b.GetName() != 'rec_kinematics_valid']
        if not rec_names or not paired.GetBranch('rec_kinematics_valid'):
            raise RuntimeError('Missing reconstructed schema: ' + str(source))
        mapping = [(name, name[4:]) for name in rec_names]
        for name in ('event_sequence', 'weight'):
            if paired.GetBranch(name) and name not in [out for _, out in mapping]:
                mapping.append((name, name))
        if not paired.GetBranch('event_sequence'):
            raise RuntimeError('Missing event_sequence: ' + str(source))
        if len({out for _, out in mapping}) != len(mapping):
            raise RuntimeError('Duplicate output branch names')
        for name in ['rec_kinematics_valid'] + [name for name, _ in mapping]:
            leaf = paired.GetLeaf(name)
            if not leaf or leaf.GetTypeName() != 'Double_t' or leaf.GetLeafCount() or leaf.GetLenStatic() != 1:
                raise RuntimeError('Expected scalar Double_t: ' + name)
        # Read the validity branch for every event; read reconstructed data only
        # for accepted rows. Never load the generated kinematic branches.
        paired.SetBranchStatus('*', 0)
        flag = array('d', [0])
        paired.SetBranchStatus('rec_kinematics_valid', 1)
        paired.SetBranchAddress('rec_kinematics_valid', flag)
        flag_branch = paired.GetBranch('rec_kinematics_valid')
        buffers = {}
        branches = []
        for name, out_name in mapping:
            buffers[out_name] = array('d', [0])
            paired.SetBranchStatus(name, 1)
            paired.SetBranchAddress(name, buffers[out_name])
            branches.append(paired.GetBranch(name))
        output_file = ROOT.TFile.Open(str(temporary), 'CREATE')
        if not output_file or output_file.IsZombie():
            raise RuntimeError('Cannot create temporary output')
        tree = ROOT.TTree('PhysicsEvents', 'Valid reconstructed MC candidates')
        for _, name in mapping:
            tree.Branch(name, buffers[name], name + '/D')
        tree.SetAutoFlush(-30000000)
        total = paired.GetEntries()
        selected = 0
        for entry in range(total):
            if flag_branch.GetEntry(entry) <= 0:
                raise RuntimeError(f'Flag read failure at {entry}')
            if flag[0] not in (0, 1):
                raise RuntimeError(f'Invalid reconstruction flag at {entry}: {flag[0]}')
            if flag[0] == 1:
                for branch in branches:
                    if branch.GetEntry(entry) <= 0:
                        raise RuntimeError(f'Branch read failure at {entry}')
                if tree.Fill() < 0:
                    raise RuntimeError('Output write failure')
                selected += 1
            if entry and entry % 10000000 == 0:
                print(f'[PROGRESS] {source.name} scanned={entry}/{total} selected={selected}', flush=True)
        if fingerprint(source) != identity:
            raise RuntimeError('Source changed during export')
        data = dict(identity, total_rows=total, selected_rows=selected,
                    branches=[name for _, name in mapping])
        if tree.Write() <= 0:
            raise RuntimeError('Tree write failure')
        ROOT.TNamed('reconstructed_export_metadata', json.dumps(data, sort_keys=True)).Write()
        output_file.Close()
        output_file = None
        inspect_output(ROOT, temporary, identity)
        # Hard-link publication is atomic and refuses to overwrite an existing file.
        os.link(temporary, dest)
        return f'[DONE] {dest} selected={selected}/{total}'
    finally:
        if output_file:
            output_file.Close()
        if input_file:
            input_file.Close()
        if temporary.exists():
            temporary.unlink()
        lock.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/volatile/clas12/thayward/dvcs_pass2_full_rebuild'))
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--scope', choices=('parts', 'merged', 'all'), default='parts',
                        help='parts avoids duplicating the same events from merged files')
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be positive')
    mc = args.root / 'mc'
    if not mc.is_dir():
        parser.error('Missing MC directory: ' + str(mc))
    files = sorted(p for p in mc.rglob('*.root') if not p.name.endswith(SUFFIX)
                   and '.partial' not in p.name and
                   (args.scope == 'all' or (args.scope == 'parts' and p.parent.name == 'parts')
                    or (args.scope == 'merged' and p.name.startswith('merged'))))
    if not files:
        parser.error('No input ROOT files found')
    for path in files:
        print('[PLAN]', path, '->', path.with_name(path.stem + SUFFIX))
    print(f'[FILES] {len(files)} execute={args.execute}', flush=True)
    if not args.execute:
        return 0
    failures = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(extract, str(path)): path for path in files}
        for future in as_completed(futures):
            try:
                print(future.result(), flush=True)
            except Exception as exc:
                failures.append(str(futures[future]))
                print('[FAILED]', futures[future], str(exc), flush=True)
    print(f'[SUMMARY] files={len(files)} failed={len(failures)}', flush=True)
    return int(bool(failures))


if __name__ == '__main__':
    raise SystemExit(main())
