#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#   ./simple_loop_filter.csh <MIN_CURRENT> <OUTPUT_DIR> <INPUT_DIR> <SUBDIR_1> [<SUBDIR_2> ...]
#
# Example:
#   ./simple_loop_filter.csh \
#       45 \
#       /scratch/thayward/output \
#       /cache/clas12/rg-a/production/decoded/6b.0.0 \
#       004003 004013 004014
#
# This will call:
#   /u/home/thayward/coatjava/bin/trigger-filter \
#       -c 45 \
#       -b 0x80000000 \
#       /cache/clas12/rg-a/production/decoded/6b.0.0/004003/*.hipo \
#       -o /scratch/thayward/output/004003_00001.hipo
# and similarly for 004013, 004014, etc.
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

# 4) Make sure the output directory exists
if ( ! -d "${OUTPUT_DIR}" ) then
    echo "Creating output directory: ${OUTPUT_DIR}"
    mkdir -p "${OUTPUT_DIR}"
endif

# 5) Loop over each subdirectory
foreach SUBDIR ($argv)
    echo "-------------------------------------------------------"
    echo "Processing subdir: $SUBDIR"

    # Build the full input wildcard
    # e.g. /cache/clas12/rg-a/production/decoded/6b.0.0/004003/*.hipo
    set INPUT_FILES = "${INPUT_DIR}/${SUBDIR}/*.hipo"

    # Build the output file name, e.g. /scratch/thayward/output/004003_00001.hipo
    set OUTPUT_FILE = "${OUTPUT_DIR}/${SUBDIR}_00001.hipo"

    # Call trigger-filter
    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        ${INPUT_FILES} \
        -o ${OUTPUT_FILE}

    echo "Done processing subdir: $SUBDIR"
    echo ""
#end