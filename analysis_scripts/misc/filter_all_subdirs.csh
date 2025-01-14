#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#
#   ./filter_all_subdirs_and_merge.csh \
#       <MERGED_FILE> \
#       <MIN_CURRENT> \
#       <OUTPUT_PARENT_DIR> \
#       <INPUT_PARENT_DIR> \
#       <SUBDIR_1> [<SUBDIR_2> ...]
#
# For example:
#
#   ./filter_all_subdirs_and_merge.csh \
#       /scratch/thayward/my_final_merged_output.hipo \
#       35 \
#       /scratch/thayward/output \
#       /cache/clas12/rg-a/production/decoded/6b.0.0 \
#       004025 004312
#
# Steps Performed:
#   1) For each subdirectory (e.g., 004025):
#      a) Gather .hipo files from <INPUT_PARENT_DIR>/<SUBDIR>
#      b) Filter them with trigger-filter -> produces <OUTPUT_PARENT_DIR>/<SUBDIR>_00001.hipo
#   2) After all subdirectories are processed, merge the per-subdirectory outputs into ONE file:
#      hipoutils -merge -o <MERGED_FILE> [list_of_subdir_outputs]
# ------------------------------------------------------------------------------

# 1) Check for enough arguments
if ( $#argv < 5 ) then
    echo "Usage: $0 MERGED_FILE MIN_CURRENT OUTPUT_PARENT INPUT_PARENT SUBDIR_1 [SUBDIR_2 ...]"
    exit 1
endif

# 2) Parse the arguments in the NEW order:
set MERGED_FILE   = $argv[1]
set MIN_CURRENT   = $argv[2]
set OUTPUT_PARENT = $argv[3]
set INPUT_PARENT  = $argv[4]

# 3) Shift away those first four, so remaining arguments are subdirectories
shift
shift
shift
shift

# 4) Prepare an array to hold the output files we create, so we can merge them later
set MERGE_LIST = ()

# 5) Loop over each subdirectory
foreach SUBDIR ($argv)
    echo "========================================================"
    echo "Processing subdirectory: $SUBDIR"

    # Build the full path to the subdirectory
    set FULL_PATH = "${INPUT_PARENT}/${SUBDIR}"

    # Gather all the .hipo files in that subdirectory
    # (Change '*.hipo' to '*' if you want to include all files.)
    set FILE_LIST = (`ls ${FULL_PATH}/*.hipo 2> /dev/null`)

    # If no .hipo files found, skip
    if ( "$FILE_LIST" == "" ) then
        echo "No .hipo files found in: ${FULL_PATH}. Skipping."
        continue
    endif

    # Construct the output filename
    set OUTPUT_FILE = "${OUTPUT_PARENT}/${SUBDIR}_00001.hipo"

    echo "-> Filtering these input files:"
    echo "$FILE_LIST"
    echo "-> Output file will be: $OUTPUT_FILE"

    # Make sure the output directory exists
    if ( ! -d "${OUTPUT_PARENT}" ) then
        echo "Creating output parent directory: ${OUTPUT_PARENT}"
        mkdir -p "${OUTPUT_PARENT}"
    endif

    # Run trigger-filter for this subdirectory
    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        $FILE_LIST \
        -o ${OUTPUT_FILE}

    # If the output was created, append it to MERGE_LIST
    if ( -f ${OUTPUT_FILE} ) then
        set MERGE_LIST = ( $MERGE_LIST $OUTPUT_FILE )
        echo "Output file created: $OUTPUT_FILE"
    else
        echo "Warning: output file was not created for $SUBDIR"
    endif

    echo "Done processing subdirectory: $SUBDIR"
    echo ""
#end

# 6) After finishing all subdirectories, merge them (if we have any outputs)
echo "========================================================"
if ( "$MERGE_LIST" == "" ) then
    echo "No output files were created, so there is nothing to merge."
    exit 0
endif

echo "Merging the following files into a single file:"
echo "$MERGE_LIST"
echo "Final merged file: $MERGED_FILE"

# Ensure the merged file's parent directory exists
set MERGED_DIR = `dirname $MERGED_FILE`
if ( ! -d $MERGED_DIR ) then
    echo "Creating directory for merged file: $MERGED_DIR"
    mkdir -p $MERGED_DIR
endif

# Perform the merge
hipoutils -merge -o ${MERGED_FILE} $MERGE_LIST

if ( -f ${MERGED_FILE} ) then
    echo "Merge completed successfully."
    echo "Final merged file: ${MERGED_FILE}"
else
    echo "Warning: Merge command did not produce ${MERGED_FILE}"
endif