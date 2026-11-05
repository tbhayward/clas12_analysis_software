#!/usr/bin/env python3
"""Run the eight proton-efficiency data/CLASDIS groups on one farm node.

Each HIPO file is an independent, resumable job.  The global worker cap is
hard-limited to 10. This launcher never compiles processing classes or the ROOT
converter. Tag and proton rows travel through named pipes into two trees in one
ROOT part file; no TSV files persist on disk.
"""

import argparse
import concurrent.futures
import hashlib
import json
import os
import pathlib
import re
import shlex
import subprocess
import threading

P = pathlib.Path
SCRIPTS = {
    "data": "processing_proton_efficiency_data.groovy",
    "mc": "processing_proton_efficiency_mc.groovy",
}


def digest_text(paths):
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="processing_scripts/proton_efficiency_campaign.json")
    parser.add_argument("--output-root")
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--max-events", type=int, default=None,
                        help="Per-HIPO event cap; 0 means unlimited")
    parser.add_argument("--pilot", action="store_true",
                        help="Run one HIPO per group, defaulting to 100000 events per HIPO")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--retry-stale", action="store_true",
                        help="Replace incomplete outputs only; valid marked outputs are retained")
    parser.add_argument("--converter-timeout", type=int, default=300,
                        help="Seconds to wait for the ROOT converter after Groovy exits")
    args = parser.parse_args()
    if not 1 <= args.workers <= 10:
        parser.error("--workers must be between 1 and the farm-node cap of 10")
    if args.max_events is not None and args.max_events < 0:
        parser.error("--max-events must be nonnegative")
    if args.converter_timeout < 1:
        parser.error("--converter-timeout must be positive")
    max_events = args.max_events if args.max_events is not None else (100000 if args.pilot else 0)

    root = P.cwd()
    config_path = P(args.config)
    config = json.loads(config_path.read_text())
    output_root = P(args.output_root or config["output_root"])
    common = root / "processing_scripts/proton_efficiency_common.groovy"
    wrappers = [root / "processing_scripts" / name for name in SCRIPTS.values()]
    converter_source = root / "processing_scripts/convert_proton_efficiency_stream.cpp"
    converter = root / "processing_scripts/convert_proton_efficiency_stream"
    jar = root / "processing_classes/dist/processing_classes.jar"
    required = [common, config_path, converter_source, *wrappers, jar]
    missing_code = [str(path) for path in required if not path.is_file()]
    if missing_code:
        raise SystemExit("Missing required campaign files:\n  " + "\n  ".join(missing_code))
    source_id = digest_text([common, config_path, converter_source, *wrappers])
    source_id += f":jar:{jar.stat().st_size}:{jar.stat().st_mtime_ns}"

    configured = {(job["period"], job["mode"]) for job in config["jobs"]}
    expected = {(period, mode) for period in ("fa18_inb", "fa18_out", "sp18_inb", "sp18_out")
                for mode in ("data", "mc")}
    if configured != expected or len(config["jobs"]) != 8:
        raise SystemExit("Campaign config must contain exactly the four data and four MC groups")

    jobs = []
    missing_inputs = []
    group_file_counts = {}
    for group in config["jobs"]:
        input_dir = P(group["input"])
        if not input_dir.is_dir():
            missing_inputs.append(str(input_dir))
            continue
        files = sorted(input_dir.rglob("*.hipo"))
        if not files:
            missing_inputs.append(str(input_dir) + " (no HIPO files)")
            continue
        if args.pilot:
            files = files[:1]
        key = (group["period"], group["mode"])
        group_file_counts[key] = len(files)
        for hipo in files:
            suffix = hashlib.sha256(str(hipo.resolve()).encode()).hexdigest()[:12]
            prefix = output_root / group["period"] / group["mode"] / "parts" / f"{hipo.stem}_{suffix}"
            jobs.append((group, hipo, prefix))

    print(f"[PLAN] groups={len(group_file_counts)}/8 files={len(jobs)} workers={args.workers} "
          f"pilot={args.pilot} max_events={max_events} execute={args.execute}", flush=True)
    for key in sorted(group_file_counts):
        print(f"[GROUP] {key[0]} {key[1]} files={group_file_counts[key]}", flush=True)
    if missing_inputs:
        print("[MISSING_INPUTS]\n  " + "\n  ".join(missing_inputs), flush=True)
    if not args.execute:
        print("[DRY_RUN] Add --execute after checking the eight groups.", flush=True)
        return
    if missing_inputs:
        raise SystemExit("Refusing to start an incomplete eight-group campaign")
    if not converter.is_file() or not os.access(converter, os.X_OK):
        raise SystemExit("Missing executable ROOT stream converter: " + str(converter) +
                         "\nCompile it before executing the campaign.")

    counts = {"done": 0, "skipped": 0, "failed": 0}
    failures = []
    lock = threading.Lock()

    def run(job):
        group, hipo, prefix = job
        prefix.parent.mkdir(parents=True, exist_ok=True)
        root_output = P(str(prefix) + ".root")
        marker_path = P(str(prefix) + ".ok.json")
        log_path = P(str(prefix) + ".processing.log")
        converter_log = P(str(prefix) + ".converter.log")
        required_marker = {
            "input": str(hipo),
            "input_size": hipo.stat().st_size,
            "input_mtime_ns": hipo.stat().st_mtime_ns,
            "period": group["period"],
            "mode": group["mode"],
            "energy": group["energy"],
            "max_events": max_events,
            "source_id": source_id,
        }
        valid = False
        if marker_path.is_file() and root_output.is_file() and root_output.stat().st_size > 0:
            try:
                marker = json.loads(marker_path.read_text())
                valid = all(marker.get(key) == value for key, value in required_marker.items())
            except (OSError, ValueError):
                valid = False
            if valid:
                checked = subprocess.run([str(converter), "--validate", str(root_output)],
                                         capture_output=True, text=True)
                rows = re.findall(r"\[PROTON_EFFICIENCY_ROOT_VALID\] tags=(\d+) protons=(\d+)",
                                  checked.stdout)
                valid = (checked.returncode == 0 and len(rows) == 1 and
                         int(rows[0][0]) == marker.get("selected_tags") and
                         int(rows[0][1]) == marker.get("proton_rows"))
        if valid:
            with lock:
                counts["done"] += 1
                counts["skipped"] += 1
            print("[SKIP_VALID]", group["period"], group["mode"], hipo.name, flush=True)
            return

        stream_prefix = P(str(prefix) + ".stream")
        tag_fifo = P(str(stream_prefix) + ".tags.tsv")
        proton_fifo = P(str(stream_prefix) + ".protons.tsv")
        temporary_root = P(str(prefix) + ".partial.root")
        legacy_tag = P(str(prefix) + ".tags.tsv")
        legacy_proton = P(str(prefix) + ".protons.tsv")
        stale = [path for path in (root_output, marker_path, tag_fifo, proton_fifo,
                                   temporary_root, legacy_tag, legacy_proton)
                 if path.exists()]
        if stale and not args.retry_stale:
            raise RuntimeError("Incomplete/stale output; inspect it or use --retry-stale: " + str(prefix))
        for path in stale:
            if path.is_dir():
                raise RuntimeError("Unexpected directory at output path: " + str(path))
            path.unlink()

        os.mkfifo(tag_fifo)
        os.mkfifo(proton_fifo)

        script = root / "processing_scripts" / SCRIPTS[group["mode"]]
        run_override = "0" if group["mode"] == "data" else "11"
        qadb_override = "0" if group["mode"] == "data" else "1"
        command_args = [
            str(root / "coatjava/bin/run-groovy"), "-cp", str(jar), str(script),
            str(hipo), str(stream_prefix), "1", str(group["energy"]),
            run_override, qadb_override, str(max_events),
        ]
        shell_command = (
            "source source_file.txt; if ($status != 0) exit 90; "
            "module load qadb/3.4.1; if ($status != 0) exit 91; exec "
            + " ".join(shlex.quote(part) for part in command_args)
        )
        converted = None
        accepted = False
        try:
            with converter_log.open("w") as conversion_log:
                converted = subprocess.Popen([str(converter), str(tag_fifo), str(proton_fifo),
                                              str(temporary_root)],
                                             stdout=conversion_log, stderr=subprocess.STDOUT)
                with log_path.open("w") as log:
                    result = subprocess.run(["/bin/tcsh", "-c", shell_command], cwd=root,
                                            stdout=log, stderr=subprocess.STDOUT)
                if result.returncode:
                    converted.terminate()
                try:
                    converter_code = converted.wait(timeout=args.converter_timeout)
                except subprocess.TimeoutExpired:
                    converted.kill()
                    converter_code = converted.wait()
            if result.returncode or converter_code:
                raise RuntimeError(f"processor={result.returncode} converter={converter_code}; "
                                   f"see {log_path} and {converter_log}")
            text = log_path.read_text(errors="replace")
            completed = re.findall(r"Complete: events=(\d+), tags=(\d+)", text)
            converted_rows = re.findall(r"\[PROTON_EFFICIENCY_ROOT_ROWS\] tags=(\d+) protons=(\d+)",
                                        converter_log.read_text(errors="replace"))
            if len(completed) != 1 or len(converted_rows) != 1:
                raise RuntimeError("Missing or ambiguous completion records for " + str(hipo))
            events, tags = map(int, completed[0])
            root_tags, root_protons = map(int, converted_rows[0])
            if tags != root_tags:
                raise RuntimeError(f"Tag-row mismatch for {hipo}: Groovy={tags}, ROOT={root_tags}")
            if max_events and events > max_events:
                raise RuntimeError(f"Event cap exceeded for {hipo}: {events}>{max_events}")
            if not temporary_root.is_file() or temporary_root.stat().st_size == 0:
                raise RuntimeError("Missing ROOT output: " + str(temporary_root))
            checked = subprocess.run([str(converter), "--validate", str(temporary_root)],
                                     capture_output=True, text=True)
            validated_rows = re.findall(r"\[PROTON_EFFICIENCY_ROOT_VALID\] tags=(\d+) protons=(\d+)",
                                        checked.stdout)
            if (checked.returncode or len(validated_rows) != 1 or
                tuple(map(int, validated_rows[0])) != (root_tags, root_protons)):
                raise RuntimeError("ROOT validation failed for " + str(temporary_root) +
                                   ": " + checked.stderr.strip())
            os.replace(temporary_root, root_output)
            marker_path.write_text(json.dumps({**required_marker, "processed_events": events,
                                               "selected_tags": tags,
                                               "proton_rows": root_protons}, sort_keys=True) + "\n")
            accepted = True
        finally:
            if converted is not None and converted.poll() is None:
                converted.kill()
                converted.wait()
            for fifo in (tag_fifo, proton_fifo):
                if fifo.exists():
                    fifo.unlink()
            if not accepted and temporary_root.exists():
                # Preserve neither a partial ROOT file nor a misleading final output.
                temporary_root.unlink()
        with lock:
            counts["done"] += 1
        print("[DONE]", group["period"], group["mode"], hipo.name,
              f"events={events}", f"tags={tags}", f"protons={root_protons}", flush=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run, job) for job in jobs]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as error:
                with lock:
                    counts["failed"] += 1
                failures.append(str(error))
                print("[FAILED]", error, flush=True)

    summary = {
        "groups": {f"{period}/{mode}": count
                   for (period, mode), count in sorted(group_file_counts.items())},
        "files": len(jobs),
        "completed": counts["done"],
        "skipped": counts["skipped"],
        "failed": counts["failed"],
        "pilot": args.pilot,
        "max_events": max_events,
        "source_id": source_id,
    }
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "campaign_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"[FINAL] completed={counts['done']}/{len(jobs)} skipped={counts['skipped']} "
          f"failed={counts['failed']}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
