#!/usr/bin/env python3
"""Fast marker inventory; never opens or scans ROOT event trees."""

import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--config", type=Path,
                        default=Path("processing_scripts/proton_efficiency_campaign.json"))
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    incomplete = 0
    source_ids = set()
    for group in config["jobs"]:
        period, mode = group["period"], group["mode"]
        parts = args.output_root / period / mode / "parts"
        markers = {}
        for marker_path in parts.glob("*.ok.json"):
            try:
                marker = json.loads(marker_path.read_text())
            except (OSError, ValueError):
                print(f"INVALID_MARKER {marker_path}")
                continue
            markers[marker_path.name[:-8]] = marker
        expected = sorted(Path(group["input"]).glob("*.hipo"))
        complete = tags = protons = 0
        for hipo in expected:
            suffix = hashlib.sha256(str(hipo.resolve()).encode()).hexdigest()[:12]
            name = f"{hipo.stem}_{suffix}"
            marker = markers.get(name)
            root_path = parts / (name + ".root")
            if (marker and root_path.is_file() and root_path.stat().st_size > 0
                    and marker.get("input") == str(hipo)
                    and marker.get("input_size") == hipo.stat().st_size
                    and marker.get("input_mtime_ns") == hipo.stat().st_mtime_ns
                    and marker.get("period") == period and marker.get("mode") == mode
                    and marker.get("max_events") == 0):
                complete += 1
                tags += marker["selected_tags"]
                protons += marker["proton_rows"]
                source_ids.add(marker.get("source_id"))
            else:
                incomplete += 1
                print(f"NEEDS_RERUN {period}/{mode} {hipo}")
        print(f"GROUP {period}/{mode} complete={complete}/{len(expected)} "
              f"tags={tags} protons={protons}", flush=True)
    print(f"TOTAL missing_or_incomplete={incomplete}")
    for source_id in sorted(source_ids):
        print(f"SOURCE_ID {source_id}")


if __name__ == "__main__":
    main()
