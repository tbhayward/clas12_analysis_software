#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#   ./filter_all_subdirs.csh <MIN_CURRENT> <OUTPUT_PARENT> <INPUT_PARENT> <SUBDIR_1> [<SUBDIR_2> ...]
#
# Example:
#   ./filter_all_subdirs.csh 45 /scratch/thayward/output /cache/clas12/rg-a/production/decoded/6b.0.0 004003 004013
#
# For each subdirectory, we gather *.hipo files, then run trigger-filter
# with -c <MIN_CURRENT> -b 0x80000000 on those files, producing
# <OUTPUT_PARENT>/<SUBDIR>_00001.hipo
# ------------------------------------------------------------------------------

# 1) Ensure enough arguments
if ( $#argv < 4 ) then
    echo "Usage: $0 MIN_CURRENT OUTPUT_PARENT INPUT_PARENT SUBDIR_1 [SUBDIR_2 ...]"
    exit 1
endif

# 2) Parse the first three arguments
set MIN_CURRENT   = $argv[1]
set OUTPUT_PARENT = $argv[2]
set INPUT_PARENT  = $argv[3]

# 3) Shift them away, leaving subdirectories in $argv
shift
shift
shift

# 4) Loop over each subdirectory
foreach SUBDIR ($argv)
    echo "========================================================"
    echo "Processing subdirectory: $SUBDIR"

    # Build the full path, e.g. /cache/clas12/rg-a/production/decoded/6b.0.0/004013
    set FULL_PATH = "${INPUT_PARENT}/${SUBDIR}"

    # Collect all *.hipo files in that subdirectory using ls
    # IMPORTANT: In csh, '2> /dev/null' doesn't work the same as in bash, so we do:
    #   ls -1 ... >& /dev/null   (redirect both stdout + stderr) and capture in backticks
    # Then we store the result in FILE_LIST as an array of filenames.
    #
    # If no files match, the array will be empty. If at least one matches, we get them all.
    #
    # The -1 ensures one filename per line, so it splits nicely into an array.
    # The >& /dev/null hides any "No match" or "No such file" errors.
    set FILE_LIST = ( `ls -1 ${FULL_PATH}/*.hipo >& /dev/null` )

    # Now let's see how many matches we got
    if ( $#FILE_LIST == 0 ) then
        echo "No *.hipo files found in: ${FULL_PATH}. Skipping."
        continue
    endif

    echo "-> Found the following *.hipo files in ${FULL_PATH}:"
    echo "$FILE_LIST"

    # Construct output file, e.g. /scratch/thayward/output/004013_00001.hipo
    set OUTPUT_FILE = "${OUTPUT_PARENT}/${SUBDIR}_00001.hipo"

    echo "-> Output file will be: $OUTPUT_FILE"

    # Make sure the output directory exists
    if ( ! -d "${OUTPUT_PARENT}" ) then
        echo "Creating directory: ${OUTPUT_PARENT}"
        mkdir -p "${OUTPUT_PARENT}"
    endif

    # Call trigger-filter
    /u/home/thayward/coatjava/bin/trigger-filter \
        -c ${MIN_CURRENT} \
        -b 0x80000000 \
        $FILE_LIST \
        -o ${OUTPUT_FILE}

    echo "Done processing subdirectory: $SUBDIR"
    echo ""
#end