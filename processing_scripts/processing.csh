#!/bin/csh

# Set the first argument to be processing single hadrons if not provided
if ( $#argv < 1 ) then
    set arg1 = "processing_scripts/processing_single_hadrons.groovy"
    echo "Warning: First argument not provided. Using default: $arg2, which processes semi-inclusive single hadron events."
else
    set arg2 = "$2"
endif

# Set the second argument to default to the RGA Fall2018 inbending nSidis skim if not provided
if ( $#argv < 2 ) then
    set arg2 = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/"
    echo "Warning: Second argument not provided. Using default: $arg2, the RGAFa18 Inbending nSidis skim."
else
    set arg2 = "$2"
endif

cd clasqaDB/
source env.csh
cd ..
coatjava/bin/run-groovy "$arg1" "$arg2" "$3" "$4" "$5" "$6"
