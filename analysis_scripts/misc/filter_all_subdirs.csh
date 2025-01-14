#!/bin/tcsh

# ------------------------------------------------------------------------------
# Usage:
#   ./filter_all_subdirs.csh <MIN_CURRENT> <OUTPUT_FILE> <PARENT_DIR> <SUBDIR_1> [<SUBDIR_2> ...]
#
# Example:
#   ./filter_all_subdirs.csh 35 /scratch/thayward/my_output.hipo \
#       /cache/clas12/rg-a/production/decoded/6b.0.0/ 004025 004312
#
# This will run:
#   /u/home/thayward/coatjava/bin/trigger-filter \
#       -c 35 \
#       -b 0x80000000 \
#       /cache/clas12/rg-a/production/decoded/6b.0.0/004025/* \
#       /cache/clas12/rg-a/production/decoded/6b.0.0/004312/* \
#       -o /scratch/thayward/my_output.hipo
# ------------------------------------------------------------------------------

# 1) Check for minimum number of arguments
if ( $#argv < 4 ) then
    echo "Usage: $0 MIN_CURRENT OUTPUT_FILE PARENT_DIR SUBDIR_1 [SUBDIR_2 ...]"
    exit 1
endif

# 2) Parse arguments
set MIN_CURRENT  = $argv[1]
set OUTPUT_FILE  = $argv[2]
set PARENT_DIR   = $argv[3]

# 3) Shift so that remaining arguments are the subdirectories
shift
shift
shift

# 4) Combine all the subdirectories' files into a single list
set FILE_LIST = ()
foreach SUBDIR ($argv)
    set FULL_PATH = "${PARENT_DIR}/${SUBDIR}"

    # Collect all files from each subdirectory
    # The wildcard * will match everything, so you may want *.hipo if you only want .hipo files.
    # But the user specifically said "the code will handle the subdirectories inside them," 
    # so let's gather everything. Adjust to *.hipo if you only want .hipo files.
    set SUBDIR_FILES = (`ls ${FULL_PATH}/*`)
    if ( "$SUBDIR_FILES" == "" ) then
        echo "Warning: no files found in ${FULL_PATH}/*"
    else
        # Append these files to FILE_LIST
        set FILE_LIST = ( $FILE_LIST $SUBDIR_FILES )
    endif
end

# 5) Check if we have at least one file
if ( "$FILE_LIST" == "" ) then
    echo "No files collected from subdirectories. Exiting."
    exit 1
endif

# 6) Call trigger-filter with the entire set of files
echo "Calling trigger-filter with the following files:"
echo "$FILE_LIST"
echo "Output: $OUTPUT_FILE"

# If you have some other arguments like -n <number> to pass, 
# add them similarly in the command below.
# E.g., -n 20011 if desired.
/u/home/thayward/coatjava/bin/trigger-filter \
    -c ${MIN_CURRENT} \
    -b 0x80000000 \
    $FILE_LIST \
    -o ${OUTPUT_FILE}