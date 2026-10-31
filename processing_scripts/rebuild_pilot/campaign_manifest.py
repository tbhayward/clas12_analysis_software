"""Frozen campaign inventory and completion verification (no ROOT mutation)."""
import hashlib
import json
import pathlib

PERIODS = ('sp18_inb', 'sp18_out', 'fa18_inb', 'fa18_out', 'sp19_inb')
MODES = ('dvcs_data', 'pi0_data', 'dvcs_rec', 'pi0_rec')

def inventory(config, cap):
    records = []
    seen = set()
    coverage = set()
    for job in config['jobs']:
        if job['period'] not in PERIODS or job['mode'] not in MODES:
            continue
        coverage.add((job['period'], job['mode']))
        folder = pathlib.Path(job['input'])
        files = sorted(folder.rglob('*.hipo')) if folder.is_dir() else []
        if not files:
            raise ValueError('Missing or empty input: ' + str(folder))
        modes = {'dvcs_rec': ('dvcs_mc',), 'pi0_rec': ('pi0_mc', 'pi0_as_dvcs')}.get(job['mode'], (job['mode'],))
        for mode in modes:
            for source in files:
                path = str(source.resolve())
                key = (mode, path)
                if key in seen:
                    raise ValueError('Duplicate input assignment: ' + str(key))
                seen.add(key)
                records.append(dict(period=job['period'], mode=mode, input=path,
                                    input_size=source.stat().st_size,
                                    input_mtime_ns=source.stat().st_mtime_ns,
                                    energy=job['energy'], current=job.get('current', 'nominal')))
    missing = set((p, m) for p in PERIODS for m in MODES) - coverage
    if missing:
        raise ValueError('Missing required period/mode groups: ' + str(sorted(missing)))
    return dict(version=1, max_events=cap, records=sorted(records, key=lambda r: (r['period'], r['mode'], r['input'])))

def freeze(base, manifest):
    payload = json.dumps(manifest, sort_keys=True, indent=2) + '\n'
    path = base / 'campaign_manifest.json'
    try:
        with path.open('x') as stream:
            stream.write(payload)
    except FileExistsError:
        if path.read_text() != payload:
            raise ValueError('Frozen campaign manifest differs; use a new output root')
    return hashlib.sha256(payload.encode()).hexdigest()

def verify(base, manifest, node=None):
    expected = []
    for record in manifest['records']:
        owner = int.from_bytes(hashlib.sha256((record['mode']+'\0'+record['input']).encode()).digest()[:8], 'big') % 2 + 1
        if node is not None and owner != node:
            continue
        expected.append(record)
    found = {}
    for marker in base.glob('**/parts/*.root.ok'):
        data = json.loads(marker.read_text())
        key = (data.get('mode'), str(pathlib.Path(data['input']).resolve()))
        if key in found:
            raise ValueError('Duplicate completion marker: ' + str(marker))
        found[key] = (marker, data)
    for record in expected:
        key = (record['mode'], record['input'])
        if key not in found:
            raise ValueError('Missing completed job: ' + str(key))
        marker, data = found[key]
        output = pathlib.Path(str(marker)[:-3])
        if not output.is_file() or not output.stat().st_size:
            raise ValueError('Missing output: ' + str(output))
        if data.get('period') != record['period'] or data.get('max_events') != manifest['max_events']:
            raise ValueError('Completion metadata mismatch: ' + str(marker))
        if data.get('input_size') != record['input_size']:
            raise ValueError('Input size mismatch: ' + str(marker))
    if node is None and set(found) != {(r['mode'], r['input']) for r in expected}:
        raise ValueError('Unexpected completed jobs outside frozen manifest')
