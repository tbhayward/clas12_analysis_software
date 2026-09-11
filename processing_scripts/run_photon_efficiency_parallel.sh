#!/usr/bin/env bash
set -euo pipefail

# Called by processing.csh for process_photon_efficiency.groovy.
# Usage:
#   run_photon_efficiency_parallel.sh INPUT OUTPUT_DIR NFILES BEAM RUN_OVERRIDE QADB_OVERRIDE IS_MC NWORKERS MX2_MIN MX2_MAX KEEP_TXT

INPUT=${1:?input HIPO file/directory required}
OUTDIR=${2:?output directory required}
NFILES=${3:-0}
BEAM=${4:-10.6041}
RUN_OVERRIDE=${5:-0}
QADB_OVERRIDE=${6:-0}
IS_MC=${7:-0}
NWORKERS=${8:-1}
MX2_MIN=${9:--1.0}
MX2_MAX=${10:-2.0}
KEEP_TXT=${11:-0}

SCRIPT="processing_scripts/process_photon_efficiency.groovy"
JAR="processing_classes/dist/processing_classes.jar"
CONVERTER_SRC="processing_scripts/convert_photon_efficiency_txt_to_root.cpp"
CONVERTER="processing_scripts/convert_photon_efficiency_txt_to_root"
mkdir -p "$OUTDIR"

if ! command -v root-config >/dev/null 2>&1; then
  echo "ERROR: root-config not found. Load/source ROOT before running." >&2
  exit 2
fi

g++ -O2 $(root-config --cflags) "$CONVERTER_SRC" -o "$CONVERTER" $(root-config --libs)

listfile=$(mktemp)
trap 'rm -f "$listfile"' EXIT
if [[ -f "$INPUT" ]]; then
  printf '%s\n' "$INPUT" > "$listfile"
else
  find "$INPUT" -type f -name '*.hipo' | sort > "$listfile"
fi
if [[ "$NFILES" =~ ^[0-9]+$ ]] && (( NFILES > 0 )); then
  head -n "$NFILES" "$listfile" > "${listfile}.limited"
  mv "${listfile}.limited" "$listfile"
fi

TOTAL=$(wc -l < "$listfile" | tr -d ' ')
if (( TOTAL == 0 )); then
  echo "ERROR: no HIPO files found under $INPUT" >&2
  exit 3
fi

echo "Photon-efficiency processing: $TOTAL HIPO files, $NWORKERS worker(s)"
echo "Output directory: $OUTDIR"
echo "Loose Mx2(ep) window: [$MX2_MIN, $MX2_MAX] GeV^2"

export OUTDIR BEAM RUN_OVERRIDE QADB_OVERRIDE IS_MC MX2_MIN MX2_MAX KEEP_TXT SCRIPT JAR CONVERTER

worker='\
hipo="$1"; \
base=$(basename "$hipo" .hipo); \
txt="$OUTDIR/${base}_photon_efficiency.txt"; \
root="$OUTDIR/${base}_photon_efficiency.root"; \
echo "[START] $hipo"; \
coatjava/bin/run-groovy -cp "$JAR" "$SCRIPT" "$hipo" "$txt" "$BEAM" "$RUN_OVERRIDE" "$QADB_OVERRIDE" "$IS_MC" "$MX2_MIN" "$MX2_MAX" && \
"$CONVERTER" "$txt" "$root" && \
{ if [[ "$KEEP_TXT" == "0" ]]; then rm -f "$txt"; fi; } && \
echo "[DONE ] $root"\
'

# -n1 gives one HIPO file to each worker.  -P supplies immediate parallelism.
xargs -d '\n' -n 1 -P "$NWORKERS" bash -c "$worker" _ < "$listfile"

echo "All requested photon-efficiency files completed."
