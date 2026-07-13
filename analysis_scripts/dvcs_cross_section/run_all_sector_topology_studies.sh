#!/usr/bin/env bash

set -euo pipefail


# ============================================================================
# User-configurable settings
# ============================================================================

MAIN_CPP="main.cpp"

# Command used to rebuild after changing main.cpp.
BUILD_CMD=(make -j8)

# Command used to run the complete DVCS analysis.
RUN_CMD=(./main)

# Final CSV produced by each analysis run.
FINAL_CSV="output/csvs/dvcs_pass2_analysis.csv"

# Permanent destination for the selected CSV from each run.
DEST_DIR="/u/home/thayward/dvcs_analysis_csvs"

# Directory for per-configuration terminal logs.
LOG_DIR="${DEST_DIR}/logs"


# ============================================================================
# Initial setup
# ============================================================================

if [[ ! -f "${MAIN_CPP}" ]]; then
    echo "ERROR: Could not find ${MAIN_CPP}."
    echo "Run this script from the directory containing main.cpp."
    exit 1
fi

mkdir -p "${DEST_DIR}"
mkdir -p "${LOG_DIR}"

BACKUP_MAIN_CPP="$(mktemp "${TMPDIR:-/tmp}/main.cpp.backup.XXXXXX")"
cp "${MAIN_CPP}" "${BACKUP_MAIN_CPP}"


# ============================================================================
# Restore main.cpp when the script exits
# ============================================================================

restore_main_cpp() {
    local exit_status=$?

    echo
    echo "Restoring original ${MAIN_CPP}..."

    if [[ -f "${BACKUP_MAIN_CPP}" ]]; then
        cp "${BACKUP_MAIN_CPP}" "${MAIN_CPP}"
        rm -f "${BACKUP_MAIN_CPP}"
    fi

    if [[ ${exit_status} -eq 0 ]]; then
        echo "Completed all requested configurations."
        echo "Original ${MAIN_CPP} has been restored."
    else
        echo "The script exited with status ${exit_status}."
        echo "Original ${MAIN_CPP} has still been restored."
    fi

    exit "${exit_status}"
}

trap restore_main_cpp EXIT
trap 'exit 130' INT
trap 'exit 143' TERM


# ============================================================================
# Function for changing the study configuration in main.cpp
# ============================================================================

set_configuration() {
    local enable_topology_filter="$1"
    local required_detector1="$2"
    local required_detector2="$3"

    local enable_electron_fd_sector_filter="$4"
    local electron_fd_sector="$5"

    local enable_proton_fd_sector_filter="$6"
    local proton_fd_sector="$7"

    local enable_proton_cd_sector_filter="$8"
    local proton_cd_sector="$9"

    local enable_photon_fd_sector_filter="${10}"
    local photon_fd_sector="${11}"

    python3 - \
        "${MAIN_CPP}" \
        "${enable_topology_filter}" \
        "${required_detector1}" \
        "${required_detector2}" \
        "${enable_electron_fd_sector_filter}" \
        "${electron_fd_sector}" \
        "${enable_proton_fd_sector_filter}" \
        "${proton_fd_sector}" \
        "${enable_proton_cd_sector_filter}" \
        "${proton_cd_sector}" \
        "${enable_photon_fd_sector_filter}" \
        "${photon_fd_sector}" <<'PYTHON'
import pathlib
import re
import sys


main_cpp = pathlib.Path(sys.argv[1])

values = {
    "global_cfg.enable_topology_filter": sys.argv[2],
    "global_cfg.required_detector1": sys.argv[3],
    "global_cfg.required_detector2": sys.argv[4],

    "global_cfg.enable_electron_fd_sector_filter": sys.argv[5],
    "global_cfg.electron_fd_sector": sys.argv[6],

    "global_cfg.enable_proton_fd_sector_filter": sys.argv[7],
    "global_cfg.proton_fd_sector": sys.argv[8],

    "global_cfg.enable_proton_cd_sector_filter": sys.argv[9],
    "global_cfg.proton_cd_sector": sys.argv[10],

    "global_cfg.enable_photon_fd_sector_filter": sys.argv[11],
    "global_cfg.photon_fd_sector": sys.argv[12],
}

text = main_cpp.read_text()

for variable, value in values.items():
    pattern = re.compile(
        rf"^(\s*{re.escape(variable)}\s*=\s*)"
        rf"(?:true|false|-?\d+)"
        rf"(\s*;.*)$",
        re.MULTILINE,
    )

    text, replacement_count = pattern.subn(
        lambda match: f"{match.group(1)}{value}{match.group(2)}",
        text,
        count=1,
    )

    if replacement_count != 1:
        raise RuntimeError(
            f"Expected exactly one assignment for {variable}, "
            f"but found {replacement_count}."
        )
    #endif
#endfor

main_cpp.write_text(text)
PYTHON
}


# ============================================================================
# Function for building, running and saving one configuration
# ============================================================================

run_configuration() {
    local output_name="$1"

    local enable_topology_filter="$2"
    local required_detector1="$3"
    local required_detector2="$4"

    local enable_electron_fd_sector_filter="$5"
    local electron_fd_sector="$6"

    local enable_proton_fd_sector_filter="$7"
    local proton_fd_sector="$8"

    local enable_proton_cd_sector_filter="$9"
    local proton_cd_sector="${10}"

    local enable_photon_fd_sector_filter="${11}"
    local photon_fd_sector="${12}"

    local output_csv="${DEST_DIR}/${output_name}.csv"
    local log_file="${LOG_DIR}/${output_name}.log"

    echo
    echo "============================================================================"
    echo "Starting configuration: ${output_name}"
    echo "Output CSV:             ${output_csv}"
    echo "Log file:               ${log_file}"
    echo "============================================================================"

    set_configuration \
        "${enable_topology_filter}" \
        "${required_detector1}" \
        "${required_detector2}" \
        "${enable_electron_fd_sector_filter}" \
        "${electron_fd_sector}" \
        "${enable_proton_fd_sector_filter}" \
        "${proton_fd_sector}" \
        "${enable_proton_cd_sector_filter}" \
        "${proton_cd_sector}" \
        "${enable_photon_fd_sector_filter}" \
        "${photon_fd_sector}"

    rm -f "${FINAL_CSV}"

    {
        echo "Configuration: ${output_name}"
        echo "Started:       $(date --iso-8601=seconds)"
        echo
        echo "Building with:"
        printf '  %q' "${BUILD_CMD[@]}"
        echo
        echo

        "${BUILD_CMD[@]}"

        echo
        echo "Running with:"
        printf '  %q' "${RUN_CMD[@]}"
        echo
        echo

        "${RUN_CMD[@]}"

        echo
        echo "Finished: $(date --iso-8601=seconds)"
    } 2>&1 | tee "${log_file}"

    if [[ ! -f "${FINAL_CSV}" ]]; then
        echo "ERROR: The expected CSV was not produced:"
        echo "       ${FINAL_CSV}"
        exit 1
    fi

    cp "${FINAL_CSV}" "${output_csv}"

    if [[ ! -s "${output_csv}" ]]; then
        echo "ERROR: Copied CSV is missing or empty:"
        echo "       ${output_csv}"
        exit 1
    fi

    echo
    echo "Saved ${output_name}:"
    echo "  ${output_csv}"
}


# ============================================================================
# Run all configurations
#
# Argument order after the output name:
#
#   topology toggle
#   required_detector1
#   required_detector2
#   electron-sector toggle
#   electron sector
#   proton-FD-sector toggle
#   proton FD sector
#   proton-CD-sector toggle
#   proton CD sector
#   photon-FD-sector toggle
#   photon FD sector
# ============================================================================


# ----------------------------------------------------------------------------
# 1. All data
# ----------------------------------------------------------------------------

run_configuration \
    "all_data" \
    false 2 0 \
    false 1 \
    false 1 \
    false 1 \
    false 1


# ----------------------------------------------------------------------------
# 2. Individual topology studies
#
# CD-FT: CD proton, FT photon -> detector1 = 2, detector2 = 0
# CD-FD: CD proton, FD photon -> detector1 = 2, detector2 = 1
# FD-FD: FD proton, FD photon -> detector1 = 1, detector2 = 1
# ----------------------------------------------------------------------------

run_configuration \
    "CD-FD" \
    true 2 1 \
    false 1 \
    false 1 \
    false 1 \
    false 1

run_configuration \
    "CD-FT" \
    true 2 0 \
    false 1 \
    false 1 \
    false 1 \
    false 1

run_configuration \
    "FD-FD" \
    true 1 1 \
    false 1 \
    false 1 \
    false 1 \
    false 1


# ----------------------------------------------------------------------------
# 3. Electron FD sectors 1 through 6
# ----------------------------------------------------------------------------

for sector in {1..6}; do
    run_configuration \
        "elec_sec${sector}" \
        false 2 0 \
        true "${sector}" \
        false 1 \
        false 1 \
        false 1
done
#endfor


# ----------------------------------------------------------------------------
# 4. Photon FD sectors 1 through 6
# ----------------------------------------------------------------------------

for sector in {1..6}; do
    run_configuration \
        "phot_sec${sector}" \
        false 2 0 \
        false 1 \
        false 1 \
        false 1 \
        true "${sector}"
done
#endfor


# ----------------------------------------------------------------------------
# 5. Proton FD sectors 1 through 6
# ----------------------------------------------------------------------------

for sector in {1..6}; do
    run_configuration \
        "prot_sec${sector}" \
        false 2 0 \
        false 1 \
        true "${sector}" \
        false 1 \
        false 1
done
#endfor


# ----------------------------------------------------------------------------
# 6. Proton CD sectors 1 through 3
# ----------------------------------------------------------------------------

for sector in {1..3}; do
    run_configuration \
        "prot_CD_sec${sector}" \
        false 2 0 \
        false 1 \
        false 1 \
        true "${sector}" \
        false 1
done
#endfor


echo
echo "============================================================================"
echo "All 25 configurations completed successfully."
echo "CSV files are stored in:"
echo "  ${DEST_DIR}"
echo
echo "Logs are stored in:"
echo "  ${LOG_DIR}"
echo "============================================================================"