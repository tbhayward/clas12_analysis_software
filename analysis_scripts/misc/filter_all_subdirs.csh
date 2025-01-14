#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#   ./filter_all_subdirs.csh <MERGED_FILE> <MIN_CURRENT> <OUTPUT_PARENT> <INPUT_PARENT> <DIR_1> [DIR_2 ...]
#
# Example:
#   ./filter_all_subdirs.csh \
#       /scratch/thayward/final_merged_output.hipo \
#       35 \
#       /scratch/thayward \
#       /cache/clas12/rg-a/production/decoded/6b.0.0 \
#       004025 004312
#
# This script:
#   1) Loops over each specified subdirectory
#       a) Finds all *.hipo files
#       b) Calls trigger-filter, writing one output per subdirectory
#   2) Merges all subdirectory outputs (OUTPUT_PARENT/*.hipo) into one file:
#       hipoutils -merge -o <MERGED_FILE> OUTPUT_PARENT/*.hipo
#   3) Removes those intermediate subdirectory output files from OUTPUT_PARENT
# ------------------------------------------------------------------------------
 
# Check for minimum number of arguments
# We need at least 5: MERGED_FILE, MIN_CURRENT, OUTPUT_PARENT, INPUT_PARENT, and at least 1 directory
if ( $#argv < 5 ) then
    echo "Usage: $0 MERGED_FILE MIN_CURRENT OUTPUT_PARENT INPUT_PARENT DIR_1 [DIR_2 ...]"
    exit 1
endif

# 1) The merged output file (new first argument)
set MERGED_FILE   = $argv[1]

# 2) The minimum current (goes to -c)
set MIN_CURRENT   = $argv[2]

# 3) The output base directory (where subdirectory outputs go)
set OUTPUT_PARENT = $argv[3]

# 4) The input base directory (where all numbered subdirs live)
set INPUT_PARENT  = $argv[4]

# Shift so $argv now holds only the subdirectories
shift
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
    set FILE_LIST = (`ls ${FULL_PATH}/*.hipo 2> /dev/null`)

    # If no files found, skip
    if ( "$FILE_LIST" == "" ) then
        echo "No *.hipo files found in ${FULL_PATH}, skipping."
        continue
    endif

    # Construct the output filename, e.g.: /scratch/thayward/004025_00001.hipo
    set OUTPUT_FILE = "${OUTPUT_PARENT}/${SUBDIR}_00001.hipo"

    # Run the filter
    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        $FILE_LIST \
        -o ${OUTPUT_FILE}

    echo "Output written to: ${OUTPUT_FILE}"
end

# After all subdirectories are processed, merge their outputs:
echo "------------------------------"
echo "Merging all subdir output files in ${OUTPUT_PARENT} into:"
echo "  ${MERGED_FILE}"
hipoutils -merge -o "${MERGED_FILE}" ${OUTPUT_PARENT}/*.hipo

# Finally, remove the intermediate subdirectory output files
echo "Removing intermediate files from: ${OUTPUT_PARENT}/*.hipo"
rm ${OUTPUT_PARENT}/*.hipo
echo "All done!"