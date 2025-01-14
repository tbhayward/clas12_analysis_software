#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#
#   ./filter_all_subdirs.csh <MIN_CURRENT> <OUTPUT_PARENT_DIR> <INPUT_PARENT_DIR> <SUBDIR_1> [<SUBDIR_2> ...]
#
# Example:
#   ./filter_all_subdirs.csh 35 /scratch/thayward /cache/clas12/rg-a/production/decoded/6b.0.0 004025 004312
#
# This produces one output file per subdirectory. For example:
#   /scratch/thayward/004025_00001.hipo
#   /scratch/thayward/004312_00001.hipo
# ------------------------------------------------------------------------------

# 1) Check for enough arguments
if ( $#argv < 4 ) then
    echo "Usage: $0 MIN_CURRENT OUTPUT_PARENT INPUT_PARENT DIR_1 [DIR_2 ...]"
    exit 1
endif

# 2) Parse arguments
set MIN_CURRENT   = $argv[1]
set OUTPUT_PARENT = $argv[2]
set INPUT_PARENT  = $argv[3]

# 3) Shift them so remaining arguments are subdirectories
shift
shift
shift

# 4) For each subdirectory:
foreach SUBDIR ($argv)
    echo "========================================================"
    echo "Processing subdirectory: $SUBDIR"

    # Full path to the subdirectory we want to read from:
    set FULL_PATH = "${INPUT_PARENT}/${SUBDIR}"

    # Check if there is at least one matching *.hipo file
    if ( -f ${FULL_PATH}/*.hipo ) then
        # Gather those files into FILE_LIST
        set FILE_LIST = ( ${FULL_PATH}/*.hipo )
    else
        echo "No .hipo files found in: ${FULL_PATH}. Skipping."
        continue
    endif

    # Construct the output file name, e.g.: /scratch/thayward/004025_00001.hipo
    set OUTPUT_FILE = "${OUTPUT_PARENT}/${SUBDIR}_00001.hipo"

    echo "-> Filtering these input files:"
    echo "$FILE_LIST"
    echo "-> Output file will be: $OUTPUT_FILE"

    # Make sure the output directory exists
    if ( ! -d "${OUTPUT_PARENT}" ) then
        echo "Creating output directory: ${OUTPUT_PARENT}"
        mkdir -p "${OUTPUT_PARENT}"
    endif

    # Run trigger-filter on the collected files
    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        $FILE_LIST \
        -o ${OUTPUT_FILE}

    echo "Done processing subdirectory: $SUBDIR"
    echo ""
#endfor