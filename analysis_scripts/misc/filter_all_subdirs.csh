#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#   ./filter_all_subdirs.csh <MERGED_FILE> <MIN_CURRENT> <OUTPUT_PARENT> <INPUT_PARENT> <SUBDIR_1> [SUBDIR_2 ...]
#
# Example:
#   ./filter_all_subdirs.csh \
#       /scratch/thayward/test.hipo \
#       35 \
#       /scratch/thayward/output \
#       /cache/clas12/rg-a/production/decoded/6b.0.0 \
#       004025 004312
#
# Steps:
#   1) For each subdirectory, gather *.hipo and call trigger-filter -> produce <SUBDIR>_00001.hipo
#   2) Merge all <OUTPUT_PARENT>/*.hipo into one file <MERGED_FILE>
#   3) Remove the intermediate <OUTPUT_PARENT>/*.hipo files
# ------------------------------------------------------------------------------

# 1) Require at least 5 arguments
if ( $#argv < 5 ) then
    echo "Usage: $0 MERGED_FILE MIN_CURRENT OUTPUT_PARENT INPUT_PARENT DIR_1 [DIR_2 ...]"
    exit 1
endif

# 2) Parse the arguments
set MERGED_FILE   = $argv[1]
set MIN_CURRENT   = $argv[2]
set OUTPUT_PARENT = $argv[3]
set INPUT_PARENT  = $argv[4]

# 3) Shift so remaining arguments are subdirectories
shift
shift
shift
shift

# 4) Loop over each subdirectory
foreach SUBDIR ($argv)
    echo "-------------------------------------------------------"
    echo "Processing directory: $SUBDIR"

    # Build full path, e.g.: /cache/clas12/rg-a/production/decoded/6b.0.0/004025
    set FULL_PATH = "${INPUT_PARENT}/${SUBDIR}"

    # Gather *.hipo with 'ls -1' (one filename per line)
    # redirect both stdout+stderr in csh with '>& /dev/null'
    set FILE_LIST = ( `ls -1 ${FULL_PATH}/*.hipo >& /dev/null` )

    # If FILE_LIST is empty, skip
    if ( $#FILE_LIST == 0 ) then
        echo "No *.hipo files found in ${FULL_PATH}, skipping."
        continue
    endif

    echo "Found files:"
    echo "$FILE_LIST"

    # Construct output file, e.g.: /scratch/thayward/output/004025_00001.hipo
    set OUTPUT_FILE = "${OUTPUT_PARENT}/${SUBDIR}_00001.hipo"

    # Ensure output parent directory exists
    if ( ! -d "${OUTPUT_PARENT}" ) then
        echo "Creating directory: ${OUTPUT_PARENT}"
        mkdir -p "${OUTPUT_PARENT}"
    endif

    # Run trigger-filter
    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        $FILE_LIST \
        -o ${OUTPUT_FILE}

    echo "Wrote filtered file to: ${OUTPUT_FILE}"
    echo ""
end

# 5) Merge step:
#    Only do so if we actually have *.hipo files in OUTPUT_PARENT
#    (Otherwise "No match" errors can occur.)
echo "-------------------------------------------------------"
echo "Merging all subdirectory outputs (if any) into: ${MERGED_FILE}"

# Check if we have any *.hipo in OUTPUT_PARENT
set MERGE_CANDIDATES = ( `ls -1 ${OUTPUT_PARENT}/*.hipo >& /dev/null` )
if ( $#MERGE_CANDIDATES == 0 ) then
    echo "No subdirectory outputs found in ${OUTPUT_PARENT} to merge. Exiting."
    exit 0
endif

echo "Merging these files:"
echo "$MERGE_CANDIDATES"

hipoutils -merge -o "${MERGED_FILE}" $MERGE_CANDIDATES

# 6) Remove intermediate files
echo "Removing intermediate files: ${OUTPUT_PARENT}/*.hipo"
rm ${OUTPUT_PARENT}/*.hipo

echo "Done! Final merged file is: ${MERGED_FILE}"