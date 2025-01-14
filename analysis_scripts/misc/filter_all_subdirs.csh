#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#   ./filter_all_subdirs.csh <MIN_CURRENT> <OUTPUT_PARENT_DIR> <INPUT_PARENT_DIR> <SUBDIR_1> [<SUBDIR_2> ...]
#
# Example:
#   ./filter_all_subdirs.csh 35 /scratch/thayward /cache/clas12/rg-a/production/decoded/6b.0.0 004025 004312
#
# This will run trigger-filter once per subdirectory:
#   - Gathers all *.hipo files in /cache/clas12/rg-a/production/decoded/6b.0.0/004025/
#     and writes a single output /scratch/thayward/004025_00001.hipo
#   - Then does the same for /cache/clas12/rg-a/production/decoded/6b.0.0/004312/
# ------------------------------------------------------------------------------

# 1) Check that we have enough arguments
if ( $#argv < 4 ) then
    echo "Usage: $0 MIN_CURRENT OUTPUT_PARENT INPUT_PARENT DIR_1 [DIR_2 ...]"
    exit 1
endif

# 2) Parse arguments
set MIN_CURRENT   = $argv[1]
set OUTPUT_PARENT = $argv[2]
set INPUT_PARENT  = $argv[3]

# 3) Shift so that $argv now holds only the subdirectories
shift
shift
shift

# 4) Iterate over each subdirectory
foreach SUBDIR ($argv)
    echo "========================================================"
    echo "Processing subdirectory: $SUBDIR"

    # Build the full path to the subdirectory
    set FULL_PATH = "${INPUT_PARENT}/${SUBDIR}"

    # Gather all the .hipo files in that subdirectory
    # If you prefer to pick up everything, you could use * instead of *.hipo.
    set FILE_LIST = (`ls ${FULL_PATH}/*.hipo 2> /dev/null`)

    # If no .hipo files found, skip
    if ( "$FILE_LIST" == "" ) then
        echo "No .hipo files found in: ${FULL_PATH}. Skipping."
        continue
    endif

    # Construct the output filename, e.g.  /scratch/thayward/004025_00001.hipo
    set OUTPUT_FILE = "${OUTPUT_PARENT}/${SUBDIR}_00001.hipo"

    echo "-> Filtering these files:"
    echo "$FILE_LIST"
    echo "-> Output file: $OUTPUT_FILE"

    # Run trigger-filter
    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        $FILE_LIST \
        -o ${OUTPUT_FILE}

    echo "Done processing subdirectory: $SUBDIR"
    echo ""
#endfor