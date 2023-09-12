#!/bin/csh

# Set the first argument to be processing single hadrons if not provided
if ( $#argv < 1 ) then
    set arg1 = "processing_scripts/processing_single_hadrons.groovy"
    echo "Warning: First argument not provided. Using default: $arg1, which processes semi-inclusive single hadron events."
else
    set arg1 = "$1"
endif

# Determine the third argument for ./convert_txt_to_root based on arg1
# Initialize to 0 as default value
set convert_arg3 = 0

# Set convert_arg3 based on the value of arg1
if ($arg1 == "processing_scripts/processing_single_hadrons.groovy") then
    set convert_arg3 = 1
else if ($arg1 == "processing_scripts/processing_inclusive.groovy") then
    set convert_arg3 = 0
else if ($arg1 == "processing_scripts/processing_dihadrons.groovy") then
    set convert_arg3 = 2
endif

# Set the second argument to default to the RGA Fall2018 inbending nSidis skim if not provided
if ( $#argv < 2 ) then
    set arg2 = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/"
    echo "Warning: Second argument not provided. Using default: $arg2, the RGAFa18 Inbending nSidis skim."
else
    set arg2 = "$2"
endif

cd clasqaDB/; source env.csh; cd ..;
coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3" "$4.txt" "$5" "$6"

# Run the convert_txt_to_root program
set txt_file = "$4.txt"
set root_file = "$4.root"
./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3
