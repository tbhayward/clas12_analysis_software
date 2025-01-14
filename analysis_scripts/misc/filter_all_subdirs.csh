#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage: ./filter_all_subdirs.csh min_current out_parent_dir input_parent_dir dir1 [dir2 ...]
# Example:
#   ./filter_all_subdirs.csh 35 /scratch/thayward /cache/clas12/rg-a/production/decoded/6b.0.0/ 004025 004312
#
# This script loops over the specified subdirectories and for each:
#   1) Finds all *.hipo files
#   2) Calls /u/home/thayward/coatjava/bin/trigger-filter with those files
#   3) Produces a single output file named: out_parent_dir/subdirectory_00001.hipo
# ------------------------------------------------------------------------------

# Check for minimum number of arguments
if ( $#argv < 4 ) then
    echo "Usage: $0 MIN_CURRENT OUTPUT_PARENT INPUT_PARENT DIR_1 [DIR_2 DIR_3 ...]"
    exit 1
endif

# 1) The minimum current:
set MIN_CURRENT = $argv[1]

# 2) The output base directory:
set OUTPUT_PARENT = $argv[2]

# 3) The input base directory (where all numbered subdirs live):
set INPUT_PARENT = $argv[3]

# Shift the arguments so $argv now holds only the subdirectories
shift
shift
shift

# Now, for each subdirectory:
foreach SUBDIR ($argv)
    echo "------------------------------"
    echo "Processing directory: $SUBDIR"

    # Build the full path for the subdirectory
    set FULL_PATH = "${INPUT_PARENT}/${SUBDIR}"

    # Gather all the .hipo files from that subdirectory
    set FILE_LIST = (`ls ${FULL_PATH}/*.hipo`)

    # If no files found, skip
    if ( "$FILE_LIST" == "" ) then
        echo "No *.hipo files found in ${FULL_PATH}, skipping."
        continue
    endif

    # Construct the output filename (for example: /scratch/thayward/004025_00001.hipo)
    set OUTPUT_FILE = "${OUTPUT_PARENT}/${SUBDIR}_00001.hipo"

    # Run the filter
    # Adjust this command to include any other flags you need (-n <some_number> etc.)
    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        $FILE_LIST \
        -o ${OUTPUT_FILE}

    echo "Output written to: ${OUTPUT_FILE}"
end