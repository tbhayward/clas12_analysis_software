#!/usr/bin/env bash
set -euo pipefail

# Usage: bash run_pid_exclusivity_parallel.sh PERIOD [NWORKERS]
# PERIOD = Su22 | Fa22 | Sp23
PERIOD="${1:?Usage: $0 Su22|Fa22|Sp23 [NWORKERS]}"
NWORKERS="${2:-8}"
SCRIPT="${PID_SKIM_SCRIPT:-$PWD/pid_exclusivity_skim.groovy}"
MERGER="${PID_MERGER_SCRIPT:-$PWD/merge_pid_exclusivity_workers.py}"
OUTDIR="/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration"
WORKDIR="${OUTDIR}/pid_exclusivity_workers/${PERIOD}"
mkdir -p "$WORKDIR"

case "$PERIOD" in
 Su22) INDIR="/cache/clas12/rg-c/production/summer22/pass1/10.5gev/NH3/dst/train/sidisdvcs"; OUTROOT="${OUTDIR}/rgc_su22_inb_NH3_pid_exclusivity.root" ;;
 Fa22) INDIR="/cache/clas12/rg-c/production/fall22/pass1/NH3/dst/train/sidisdvcs"; OUTROOT="${OUTDIR}/rgc_fa22_inb_NH3_pid_exclusivity.root" ;;
 Sp23) INDIR="/cache/clas12/rg-c/production/spring23/pass1/NH3/dst/train/sidisdvcs"; OUTROOT="${OUTDIR}/rgc_sp23_inb_NH3_pid_exclusivity.root" ;;
 *) echo "Unknown period: $PERIOD" >&2; exit 2 ;;
esac

# One independent JVM per HIPO. This is safer/faster than threading HipoDataSource
# inside one Groovy process and parallelizes cleanly across ifarm cores.
find "$INDIR" -maxdepth 1 -type f -name 'sidisdvcs_*.hipo' -print0 | sort -z | \
  xargs -0 -n1 -P "$NWORKERS" bash -c '
    h="$1"; base=$(basename "$h" .hipo); out="'"$WORKDIR"'/${base}.txt"
    if [[ -s "$out" ]]; then echo "SKIP existing $out"; exit 0; fi
    echo "START $h"
    run-groovy "'"$SCRIPT"'" "$h" "$out.tmp" "'"$PERIOD"'"
    mv "$out.tmp" "$out"
  ' _

python3 "$MERGER" "$WORKDIR/sidisdvcs_*.txt" "$OUTROOT"
echo "FINAL: $OUTROOT"
