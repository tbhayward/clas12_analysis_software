#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#   ./filter_all_subdirs.csh <MIN_CURRENT> <OUTPUT_PARENT> <INPUT_PARENT> <SUBDIR_1> [SUBDIR_2 ...]
#
# Example:
#   ./filter_all_subdirs.csh 45 /scratch/thayward/output /cache/clas12/rg-a/production/decoded/6b.0.0 004003 004013 ...
#
# This will:
#   1) For each <SUBDIR>:
#       - Look in <INPUT_PARENT>/<SUBDIR> for any *.hipo files
#       - If found, runs trigger-filter, writing the output to:
#           <OUTPUT_PARENT>/<SUBDIR>_00001.hipo
#   2) If no *.hipo files are found in a given subdir, the script skips it.
#
# No merging is performed in this version.
# ------------------------------------------------------------------------------

# 1) Ensure we have at least 4 arguments
if ( $#argv < 4 ) then
    echo "Usage: $0 MIN_CURRENT OUTPUT_PARENT INPUT_PARENT SUBDIR_1 [SUBDIR_2 ...]"
    exit 1
endif

# 2) Parse the first three arguments
set MIN_CURRENT   = $argv[1]
set OUTPUT_PARENT = $argv[2]
set INPUT_PARENT  = $argv[3]

# 3) Shift so that any remaining arguments are subdirectories
shift
shift
shift

# 4) Loop over each subdirectory
foreach SUBDIR ($argv)
    echo "========================================================"
    echo "Processing subdirectory: $SUBDIR"

    # Build the full path to that subdirectory
    # e.g., /cache/clas12/rg-a/production/decoded/6b.0.0/004003
    set FULL_PATH = "${INPUT_PARENT}/${SUBDIR}"

    # Check if there's at least one matching file named *.hipo
    # in that directory (using the -f test in csh).
    if ( -f ${FULL_PATH}/*.hipo ) then
        # If the test passes, we gather them:
        set FILE_LIST = ( ${FULL_PATH}/*.hipo )
    else
        echo "No .hipo files found in: ${FULL_PATH}. Skipping."
        continue
    endif

    # Construct the output file name, e.g. /scratch/thayward/output/004003_00001.hipo
    set OUTPUT_FILE = "${OUTPUT_PARENT}/${SUBDIR}_00001.hipo"

    echo "-> Filtering these files:"
    echo "$FILE_LIST"
    echo "-> Output file will be: $OUTPUT_FILE"

    # Ensure the output directory exists
    if ( ! -d "${OUTPUT_PARENT}" ) then
        echo "Creating directory: ${OUTPUT_PARENT}"
        mkdir -p "${OUTPUT_PARENT}"
    endif

    # Run the trigger-filter
    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        $FILE_LIST \
        -o ${OUTPUT_FILE}

    echo "Done processing subdirectory: $SUBDIR"
    echo ""
#end