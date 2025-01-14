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
# Example:
#   ./filter_all_subdirs_and_merge.csh \
#       /scratch/thayward/my_final_merged_output.hipo \
#       35 \
#       /scratch/thayward/output \
#       /cache/clas12/rg-a/production/decoded/6b.0.0 \
#       004025 004312
#
# Steps:
#   1) For each specified subdirectory (e.g. 004025):
#      a) Gather all .hipo files from <INPUT_PARENT_DIR>/<SUBDIR> (e.g. 004025/*.hipo)
#      b) Filter them using trigger-filter -> produce <OUTPUT_PARENT_DIR>/<SUBDIR>_00001.hipo
#   2) After all are processed, merge all such outputs into a single file:
#      hipoutils -merge -o <MERGED_FILE> [the list of per-subdir outputs]
# ------------------------------------------------------------------------------

# 1) Check we have enough arguments
if ( $#argv < 5 ) then
    echo "Usage: $0 MERGED_FILE MIN_CURRENT OUTPUT_PARENT INPUT_PARENT SUBDIR_1 [SUBDIR_2 ...]"
    exit 1
endif

# 2) Parse the arguments in the new order
set MERGED_FILE   = $argv[1]
set MIN_CURRENT   = $argv[2]
set OUTPUT_PARENT = $argv[3]
set INPUT_PARENT  = $argv[4]

# 3) Shift them out so $argv holds just the subdirectories
shift
shift
shift
shift

# 4) We'll store the newly created filtered files in MERGE_LIST
set MERGE_LIST = ()

# 5) Loop over each subdirectory
foreach SUBDIR ($argv)
    echo "========================================================"
    echo "Processing subdirectory: $SUBDIR"

    # Full path = e.g. /cache/clas12/rg-a/production/decoded/6b.0.0/004025
    set FULL_PATH = "${INPUT_PARENT}/${SUBDIR}"

    # --- In csh, "2> /dev/null" does not work. We'll either: 
    # ---  (A) redirect *both* stdout & stderr with '>& /dev/null', or 
    # ---  (B) do a simple check with '-f' or globbing.
    #
    # Here, let's just redirect both streams to /dev/null to avoid errors:
    # This means if "ls" fails, it won't print an error to the terminal.
    #
    # We'll also store the output of ls into FILE_LIST.
    # 
    # The syntax: `ls -1 ... >& /dev/null` = redirect *both* stdout and stderr in csh.
    # Then we use "echo `ls ...`" pattern to capture that listing.
    # 
    # Trick: We can do something like this:
    # set FILE_LIST = ( `ls -1 ${FULL_PATH}/*.hipo >& /dev/null` )
    #
    # If no files match, the variable will be empty.
    #
    # Alternatively, we can do a glob check:
    # if ( -f ${FULL_PATH}/*.hipo ) then
    #     set FILE_LIST = ( ${FULL_PATH}/*.hipo )
    # endif
    #
    # We'll go with the ls approach for consistency.
    #
    set FILE_LIST = ( `ls -1 ${FULL_PATH}/*.hipo >& /dev/null` )

    # If no .hipo files found, skip
    if ( "$FILE_LIST" == "" ) then
        echo "No .hipo files found in: ${FULL_PATH}. Skipping."
        continue
    endif

    # Construct the output file name, e.g.: /scratch/thayward/output/004025_00001.hipo
    set OUTPUT_FILE = "${OUTPUT_PARENT}/${SUBDIR}_00001.hipo"

    echo "-> Filtering these input files:"
    echo "$FILE_LIST"
    echo "-> Output file will be: $OUTPUT_FILE"

    # Make sure output directory exists
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

    # If the output was created, append it to our MERGE_LIST
    if ( -f ${OUTPUT_FILE} ) then
        set MERGE_LIST = ( $MERGE_LIST $OUTPUT_FILE )
        echo "Output file created: $OUTPUT_FILE"
    else
        echo "WARNING: no output file was created for $SUBDIR"
    endif

    echo "Done processing subdirectory: $SUBDIR"
    echo ""
#end

# 6) Merge everything, if we got any outputs
echo "========================================================"
if ( "$MERGE_LIST" == "" ) then
    echo "No output files were created, nothing to merge."
    exit 0
endif

echo "Merging the following files into a single file:"
echo "$MERGE_LIST"
echo "Final merged file: $MERGED_FILE"

# Make sure the merged-file directory exists
set MERGED_DIR = `dirname "${MERGED_FILE}"`
if ( ! -d "${MERGED_DIR}" ) then
    echo "Creating directory for merged file: ${MERGED_DIR}"
    mkdir -p "${MERGED_DIR}"
endif

# Merge them with hipoutils
hipoutils -merge -o "${MERGED_FILE}" $MERGE_LIST

if ( -f "${MERGED_FILE}" ) then
    echo "Merge completed successfully."
    echo "Final merged file: ${MERGED_FILE}"
else
    echo "WARNING: Merge command did not produce ${MERGED_FILE}"
endif