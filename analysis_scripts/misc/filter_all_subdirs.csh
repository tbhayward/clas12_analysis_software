#!/bin/tcsh

# ------------------------------------------------------------------------------
# Minimal Script with Debug Prints
# ------------------------------------------------------------------------------
# Usage:
#   ./filter_all_subdirs.csh <MIN_CURRENT> <OUTPUT_DIR> <INPUT_DIR> <SUBDIR_1> [<SUBDIR_2> ...]
#
# Example:
#   ./filter_all_subdirs.csh 45 /scratch/thayward/output /cache/clas12/rg-a/production/decoded/6b.0.0 004003 004013
#
# This will call trigger-filter once for each subdirectory:
#   /u/home/thayward/coatjava/bin/trigger-filter -c 45 -b 0x80000000 /path/<SUBDIR>/*.hipo -o <OUTPUT_DIR>/<SUBDIR>_00001.hipo
# ------------------------------------------------------------------------------
 
# 1) Check for enough arguments
if ( $#argv < 4 ) then
    echo "Usage: $0 MIN_CURRENT OUTPUT_DIR INPUT_DIR SUBDIR_1 [SUBDIR_2 ...]"
    exit 1
endif

# 2) Parse arguments
set MIN_CURRENT = $argv[1]
set OUTPUT_DIR  = $argv[2]
set INPUT_DIR   = $argv[3]

# 3) Shift so remaining arguments are subdirectories
shift
shift
shift

# --- DEBUG: Show how many arguments remain and what they are
echo "DEBUG: After shifting three arguments, the script sees $#argv subdirectory(ies)."
echo "DEBUG: Subdirectory list: $argv"
echo ""

# 4) Ensure the output directory exists
if ( ! -d "${OUTPUT_DIR}" ) then
    echo "Creating output directory: ${OUTPUT_DIR}"
    mkdir -p "${OUTPUT_DIR}"
endif

# 5) Loop over each subdirectory
foreach SUBDIR ($argv)
    echo "-------------------------------------------------------"
    echo "Processing subdir: $SUBDIR"

    # Build input wildcard, e.g. /cache/clas12/rg-a/production/decoded/6b.0.0/004003/*.hipo
    set INPUT_FILES = "${INPUT_DIR}/${SUBDIR}/*.hipo"
    set OUTPUT_FILE = "${OUTPUT_DIR}/${SUBDIR}_00001.hipo"

    echo "-> Running trigger-filter on: ${INPUT_FILES}"
    echo "-> Output file: ${OUTPUT_FILE}"

    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        ${INPUT_FILES} \
        -o ${OUTPUT_FILE}

    echo "Done processing subdir: $SUBDIR"
    echo ""
#end