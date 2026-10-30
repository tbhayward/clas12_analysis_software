#!/usr/bin/env python3
"""Read-only inventory of paired DVCS/AAOgen ROOT files for integration."""
from pathlib import Path
import argparse
import json

PERIODS = ('sp18_inb', 'sp18_out', 'fa18_inb', 'fa18_out', 'sp19_inb')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/volatile/clas12/thayward/dvcs_pass2_full_rebuild'))
    parser.add_argument('--output', type=Path, default=Path('rebuilt_mc_manifest.json'))
    args = parser.parse_args()
    records = []
    for period in PERIODS:
        for category in ('dvcs_mc', 'pi0_mc'):
            parent = args.root / 'mc' / period / category
            if not parent.is_dir():
                print(f'WARNING: missing {parent}')
                continue
            for current in sorted(p for p in parent.iterdir() if p.is_dir()):
                parts = current / 'parts'
                files = sorted(parts.glob('*.root')) if parts.is_dir() else []
                records.append({'period': period, 'category': category, 'configuration': current.name,
                                'directory': str(parts), 'count': len(files), 'files': [str(p) for p in files]})
                print(f'{period:10s} {category:9s} {current.name:12s} {len(files):4d} files')
    args.output.write_text(json.dumps({'root': str(args.root), 'records': records}, indent=2) + '\n')
    print(f'Wrote {args.output} ({len(records)} configurations)')

if __name__ == '__main__':
    main()
