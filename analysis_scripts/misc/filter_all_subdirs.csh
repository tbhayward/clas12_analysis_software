#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#   ./filter_all_subdirs.csh <MERGED_FILE> <MIN_CURRENT> <OUTPUT_PARENT> <INPUT_PARENT> <DIR_1> [DIR_2 ...]
#
# Example:
#   ./filter_all_subdirs.csh \
#       /scratch/thayward/test.hipo \
#       35 \
#       /scratch/thayward/output \
#       /cache/clas12/rg-a/production/decoded/6b.0.0 \
#       004025 004312
#
# What it does:
#   1) Loops over each subdirectory (e.g., 004025), finds *.hipo files, and runs trigger-filter:
#        /u/home/thayward/coatjava/bin/trigger-filter -c <MIN_CURRENT> -b 0x80000000 <files> -o <subdir_output>
#   2) After all subdirectories are processed, merges them:
#        hipoutils -merge -o <MERGED_FILE> <OUTPUT_PARENT>/*.hipo
#   3) Removes intermediate subdirectory outputs from <OUTPUT_PARENT>/*.hipo
#
# IMPORTANT:
#   This is the exact same logic as your prior "works before merging" script,
#   but now the final merged file is the *first* argument, and a merge is done at the end.
# ------------------------------------------------------------------------------

# 1) Require at least 5 arguments now:
#    merged_file, min_current, output_parent, input_parent, subdir(s)
if ( $#argv < 5 ) then
    echo "Usage: $0 MERGED_FILE MIN_CURRENT OUTPUT_PARENT INPUT_PARENT DIR_1 [DIR_2 ...]"
    exit 1
endif

# 2) Assign arguments in the new order
set MERGED_FILE   = $argv[1]
set MIN_CURRENT   = $argv[2]
set OUTPUT_PARENT = $argv[3]
set INPUT_PARENT  = $argv[4]

# 3) Shift away these four, leaving only the subdirectories in $argv
shift
shift
shift
shift

# 4) Process each subdirectory exactly as before
foreach SUBDIR ($argv)
    echo "------------------------------"
    echo "Processing directory: $SUBDIR"

    # Build the full path for the subdirectory
    set FULL_PATH = "${INPUT_PARENT}/${SUBDIR}"

    # Gather all the .hipo files from that subdirectory (as before)
    set FILE_LIST = (`ls ${FULL_PATH}/*.hipo`)

    # If no files found, skip
    if ( "$FILE_LIST" == "" ) then
        echo "No *.hipo files found in ${FULL_PATH}, skipping."
        continue
    endif

    # Construct the output filename (e.g. /scratch/thayward/004025_00001.hipo)
    set OUTPUT_FILE = "${OUTPUT_PARENT}/${SUBDIR}_00001.hipo"

    # Run the filter (same as your old version)
    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        $FILE_LIST \
        -o ${OUTPUT_FILE}

    echo "Output written to: ${OUTPUT_FILE}"
end

# 5) Now merge all of those subdirectory outputs into the first argument (MERGED_FILE)
echo "------------------------------"
echo "Merging all .hipo files from ${OUTPUT_PARENT} into: ${MERGED_FILE}"
hipoutils -merge -o "${MERGED_FILE}" ${OUTPUT_PARENT}/*.hipo

# 6) Remove the intermediate files
echo "Removing intermediate files: ${OUTPUT_PARENT}/*.hipo"
rm ${OUTPUT_PARENT}/*.hipo

echo "Done! Final merged file: ${MERGED_FILE}"