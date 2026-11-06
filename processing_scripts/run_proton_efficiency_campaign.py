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
import time

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
    parser.add_argument("--report-seconds", type=int, default=900,
                        help="Interval between percentage progress reports")
    parser.add_argument("--only-mode", choices=("data", "mc"),
                        help="Process only data or MC jobs from the complete campaign config")
    parser.add_argument("--contiguous-part",
                        help="Process one contiguous file-range part, as INDEX/TOTAL (1-based)")
    parser.add_argument("--reverse-order", action="store_true",
                        help="Run selected groups and files in reverse order")
    args = parser.parse_args()
    if not 1 <= args.workers <= 10:
        parser.error("--workers must be between 1 and the farm-node cap of 10")
    if args.max_events is not None and args.max_events < 0:
        parser.error("--max-events must be nonnegative")
    if args.converter_timeout < 1:
        parser.error("--converter-timeout must be positive")
    if args.report_seconds < 5:
        parser.error("--report-seconds must be at least 5")
    part_index, part_total = 1, 1
    if args.contiguous_part:
        match = re.fullmatch(r"(\d+)/(\d+)", args.contiguous_part)
        if not match:
            parser.error("--contiguous-part must have the form INDEX/TOTAL")
        part_index, part_total = map(int, match.groups())
        if part_total < 1 or not 1 <= part_index <= part_total:
            parser.error("--contiguous-part requires 1 <= INDEX <= TOTAL")
    max_events = args.max_events if args.max_events is not None else (100000 if args.pilot else 0)

    root = P.cwd()
    config_path = P(args.config)
    config = json.loads(config_path.read_text())
    output_root = P(args.output_root or config["output_root"])
    common = root / "processing_scripts/proton_efficiency_common.groovy"
    wrappers = [root / "processing_scripts" / name for name in SCRIPTS.values()]
    converter_source = root / "processing_scripts/convert_proton_efficiency_stream.cpp"
    converter = root / "processing_scripts/.build/convert_proton_efficiency_stream"
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
        if args.only_mode and group["mode"] != args.only_mode:
            continue
        input_dir = P(group["input"])
        if not input_dir.is_dir():
            missing_inputs.append(str(input_dir))
            continue
        files = sorted(input_dir.rglob("*.hipo"))
        if not files:
            missing_inputs.append(str(input_dir) + " (no HIPO files)")
            continue
        if args.pilot:
            # Avoid systematically sampling run-start files, whose first event
            # slices can be wholly rejected by QADB. Use a deterministic
            # interior file while retaining exactly one HIPO per group.
            files = [files[len(files) // 2]]
        elif part_total > 1:
            # Contiguous ranges let an additional node start at the far end of
            # a group while the original full-campaign node advances from the
            # beginning. Integer boundaries cover every file exactly once.
            begin = len(files) * (part_index - 1) // part_total
            end = len(files) * part_index // part_total
            files = files[begin:end]
        if not files:
            continue
        key = (group["period"], group["mode"])
        group_file_counts[key] = len(files)
        for hipo in files:
            suffix = hashlib.sha256(str(hipo.resolve()).encode()).hexdigest()[:12]
            prefix = output_root / group["period"] / group["mode"] / "parts" / f"{hipo.stem}_{suffix}"
            jobs.append((group, hipo, prefix))

    if args.reverse_order:
        jobs.reverse()
    selection = (f"mode={args.only_mode or 'all'} part={part_index}/{part_total} "
                 f"reverse={args.reverse_order}")
    print(f"[PLAN] groups={len(group_file_counts)}/8 files={len(jobs)} workers={args.workers} "
          f"pilot={args.pilot} max_events={max_events} {selection} execute={args.execute}", flush=True)
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

    total_bytes = sum(hipo.stat().st_size for _, hipo, _ in jobs)
    group_total_bytes = {}
    group_total_files = {}
    for group, hipo, _ in jobs:
        key = f"{group['period']}/{group['mode']}"
        group_total_bytes[key] = group_total_bytes.get(key, 0) + hipo.stat().st_size
        group_total_files[key] = group_total_files.get(key, 0) + 1
    counts = {"done": 0, "skipped": 0, "failed": 0, "bytes": 0,
              "events": 0, "active": {},
              "group_bytes": {key: 0 for key in group_total_bytes},
              "group_files": {key: 0 for key in group_total_files}}
    failures = []
    lock = threading.Lock()
    started = time.monotonic()

    def active_progress(item):
        processed = 0
        expected = None
        try:
            text = item["log"].read_text(errors="replace")
            progress = re.findall(r"\[HIPO_PROGRESS\].*?events_processed=(\d+)", text)
            totals = re.findall(r"number of\s+events\s*:\s*(\d+)", text)
            if progress:
                processed = int(progress[-1])
            if totals:
                expected = int(totals[-1])
                if max_events:
                    expected = min(expected, max_events)
        except OSError:
            pass
        fraction = min(1.0, processed / expected) if expected else 0.0
        return processed, expected, fraction

    stop_reporter = threading.Event()

    def reporter():
        while not stop_reporter.wait(args.report_seconds):
            with lock:
                done = counts["done"]
                skipped = counts["skipped"]
                failed = counts["failed"]
                completed_bytes = counts["bytes"]
                completed_events = counts["events"]
                active = list(counts["active"].values())
                group_bytes = dict(counts["group_bytes"])
                group_files = dict(counts["group_files"])
            active_details = []
            observed_events = completed_events
            estimated_bytes = float(completed_bytes)
            active_group_bytes = {key: 0.0 for key in group_total_bytes}
            for item in active:
                processed, expected, fraction = active_progress(item)
                observed_events += processed
                contribution = item["bytes"] * fraction
                estimated_bytes += contribution
                active_group_bytes[item["group"]] += contribution
                percent = 100.0 * fraction
                active_details.append(f"{item['label']}={percent:.1f}% ({processed}/{expected or '?'})")
            elapsed = time.monotonic() - started
            work_fraction = estimated_bytes / total_bytes if total_bytes else 1.0
            work_percent = 100.0 * min(1.0, work_fraction)
            file_percent = 100.0 * done / len(jobs) if jobs else 100.0
            rate = observed_events / elapsed if elapsed > 0 else 0.0
            eta_seconds = elapsed * (1.0 / work_fraction - 1.0) if work_fraction > 0 else None
            eta = "n/a" if eta_seconds is None else f"{eta_seconds/60.0:.1f} min"
            group_details = []
            for key in sorted(group_total_bytes):
                numerator = group_bytes[key] + active_group_bytes[key]
                percent = 100.0 * numerator / group_total_bytes[key] if group_total_bytes[key] else 100.0
                group_details.append(f"{key}={percent:.1f}%[{group_files[key]}/{group_total_files[key]}]")
            print(f"[PROGRESS] work={work_percent:.2f}% files={done}/{len(jobs)}({file_percent:.2f}%) "
                  f"active={len(active)} skipped={skipped} failed={failed} rate={rate:.1f} events/s "
                  f"elapsed={elapsed/60.0:.1f} min eta={eta}", flush=True)
            print("[GROUP_PROGRESS] " + " ".join(group_details), flush=True)
            if active_details:
                print("[ACTIVE] " + " | ".join(active_details), flush=True)

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
                counts["bytes"] += hipo.stat().st_size
                counts["events"] += int(marker.get("processed_events", 0))
                key = f"{group['period']}/{group['mode']}"
                counts["group_bytes"][key] += hipo.stat().st_size
                counts["group_files"][key] += 1
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
        group_key = f"{group['period']}/{group['mode']}"
        active_key = str(prefix)
        with lock:
            counts["active"][active_key] = {
                "label": f"{group_key}:{hipo.name}", "group": group_key,
                "bytes": hipo.stat().st_size, "log": log_path,
            }
        print("[START]", group["period"], group["mode"], hipo.name,
              f"size_gb={hipo.stat().st_size/1e9:.2f}", flush=True)
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
            with lock:
                counts["active"].pop(active_key, None)
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
            counts["bytes"] += hipo.stat().st_size
            counts["events"] += events
            counts["group_bytes"][group_key] += hipo.stat().st_size
            counts["group_files"][group_key] += 1
            completed_now = counts["done"]
        print("[DONE]", group["period"], group["mode"], hipo.name,
              f"events={events}", f"tags={tags}", f"protons={root_protons}",
              f"files={completed_now}/{len(jobs)} ({100.0*completed_now/len(jobs):.2f}%)", flush=True)

    reporter_thread = threading.Thread(target=reporter, daemon=True)
    reporter_thread.start()
    try:
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
    finally:
        stop_reporter.set()
        reporter_thread.join(timeout=2)

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
    summary["only_mode"] = args.only_mode
    summary["contiguous_part"] = f"{part_index}/{part_total}"
    summary["reverse_order"] = args.reverse_order
    if args.only_mode or part_total > 1:
        summary_name = f"campaign_summary_{args.only_mode or 'all'}_part{part_index}of{part_total}.json"
    else:
        summary_name = "campaign_summary.json"
    (output_root / summary_name).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"[FINAL] completed={counts['done']}/{len(jobs)} skipped={counts['skipped']} "
          f"failed={counts['failed']}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
