#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#   ./filter_all_subdirs.csh <MIN_CURRENT> <OUTPUT_DIR> <INPUT_DIR> <SUBDIR_1> [<SUBDIR_2> ...]
#
# Example:
#   ./filter_all_subdirs.csh \
#       45 \
#       /scratch/thayward/output \
#       /cache/clas12/rg-a/production/decoded/6b.0.0 \
#       004003 004013 004014
#
# For each subdirectory, it calls trigger-filter with:
#   -c <MIN_CURRENT>
#   -b 0x80000000
#   <INPUT_DIR>/<SUBDIR>/*.hipo
#   -o <OUTPUT_DIR>/<SUBDIR>_00001.hipo
#
# No checks, no merging—just one run per subdirectory.
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

# 4) Ensure the output directory exists
if ( ! -d "${OUTPUT_DIR}" ) then
    echo "Creating output directory: ${OUTPUT_DIR}"
    mkdir -p "${OUTPUT_DIR}"
endif

# 5) Loop over each subdirectory
foreach SUBDIR ($argv)
    echo "-------------------------------------------------------"
    echo "Processing subdir: $SUBDIR"

    # Build the wildcard for input files, e.g. /path/to/6b.0.0/004003/*.hipo
    set INPUT_FILES = "${INPUT_DIR}/${SUBDIR}/*.hipo"

    # Build the output file name, e.g. /scratch/thayward/output/004003_00001.hipo
    set OUTPUT_FILE = "${OUTPUT_DIR}/${SUBDIR}_00001.hipo"

    echo "-> Running trigger-filter on: ${INPUT_FILES}"
    echo "-> Output: ${OUTPUT_FILE}"

    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        ${INPUT_FILES} \
        -o ${OUTPUT_FILE}

    echo "Done processing subdir: $SUBDIR"
    echo ""
#end