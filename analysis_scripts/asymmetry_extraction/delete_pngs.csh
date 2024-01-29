#!/bin/csh

# Script to delete .png files in output/ and its subdirectories

echo "Are you sure you want to delete all .png files in output/ and its subdirectories? (Y/N)"
set user_input = $<

if ($user_input == "Y") then
    echo "Deleting .png files..."
    find output/ -name '*.png' -exec rm {} \;
    echo "Deletion complete."
else
    echo "Deletion aborted by user."
endif
