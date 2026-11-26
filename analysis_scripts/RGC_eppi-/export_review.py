#!/usr/bin/env python3
"""Copy compact analysis products for review, leaving event data and caches in place."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile

HERE = Path(__file__).resolve().parent
SKIP_DIRS = {'selected', 'calibration_cache', '__pycache__', 'cache', 'caches'}
SUFFIXES = {'.png', '.json', '.csv', '.txt', '.html', '.md', '.npz'}


def export_review(source, destination, run_info=None, datasets=None):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if not source.is_dir():
        raise ValueError('Analysis output directory does not exist: '+str(source))
    if source == destination or source in destination.parents or destination in source.parents:
        raise ValueError('Review destination must be separate from analysis outputs')
    destination.mkdir(parents=True, exist_ok=True)
    records = []
    for folder, dirs, files in os.walk(source, followlinks=False):
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS and not (Path(folder)/d).is_symlink())
        for name in sorted(files):
            p = Path(folder)/name
            if p.is_symlink() or p.suffix.lower() not in SUFFIXES or 'replica' in name.lower() or '.partial.' in name:
                continue
            relative = p.relative_to(source)
            target = destination/relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, target)
            records.append(dict(path=str(relative), bytes=target.stat().st_size, sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
    for name, p in [('run_info.csv', run_info or HERE/'imports/clas12_run_info_ND3_source.csv'), ('datasets.json', datasets or HERE/'imports/datasets.json')]:
        target = destination/'input_snapshots'/name
        target.parent.mkdir(exist_ok=True)
        shutil.copy2(p, target)
        records.append(dict(path=str(target.relative_to(destination)), bytes=target.stat().st_size, sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
    # Remove only files owned by an earlier export, so reruns cannot mix stale products.
    manifest = destination/'review_manifest.json'
    if manifest.exists():
        previous = json.loads(manifest.read_text())
        current = {x['path'] for x in records}
        for entry in previous.get('files', []):
            p = (destination/entry['path']).resolve()
            if destination in p.parents and entry['path'] not in current:
                p.unlink(missing_ok=True)
    manifest.write_text(json.dumps(dict(source=str(source), files=records, omitted='Selected ROOT files, calibration caches, NPY arrays and replica ensembles. Registries retain source paths and are provenance, not a portable fitting input.'), indent=2)+'\n')
    (destination/'REVIEW_README.txt').write_text('Compact review copy: figures, fit JSONs, tables, covariance/histogram NPZs, and input metadata.\nLarge event files, calibration caches and replica ensembles remain at '+str(source)+'.\nThe full source output is required to rerun fits.\n')
    archive = destination.with_name(destination.name+'.tar.gz')
    temporary = archive.with_name(archive.name+'.tmp')
    with tarfile.open(temporary, 'w:gz') as tar:
        tar.add(destination, arcname=destination.name)
    temporary.replace(archive)
    print(f'Review copy: {destination} ({sum(x["bytes"] for x in records)/1024**2:.1f} MiB)')
    print(f'Transfer archive: {archive} ({archive.stat().st_size/1024**2:.1f} MiB)', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--destination', type=Path, default=HERE/'review_loose_20261124')
    p.add_argument('--run-info', type=Path)
    p.add_argument('--datasets', type=Path)
    a = p.parse_args()
    export_review(a.source, a.destination, a.run_info, a.datasets)


if __name__ == '__main__':
    main()
