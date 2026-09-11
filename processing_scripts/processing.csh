#!/bin/csh

source source_file.txt


# Dedicated photon-efficiency path.  This mode processes one HIPO file per worker,
# converts each completed text checkpoint immediately to a ROOT tree, and can run
# multiple files in parallel.  It exits here so the legacy argument conventions
# below remain unchanged for all existing processing scripts.
if ( $#argv >= 1 && "$1" == "processing_scripts/process_photon_efficiency.groovy" ) then
    if ( $#argv < 3 ) then
        echo "Usage: ./processing_scripts/processing.csh processing_scripts/process_photon_efficiency.groovy INPUT_HIPO_OR_DIR OUTPUT_DIR [NFILES=0] [BEAM=10.6041] [RUN_OVERRIDE=0] [QADB_OVERRIDE=0] [IS_MC=0] [NWORKERS=1] [MX2_MIN=-1.0] [MX2_MAX=2.0] [KEEP_TXT=0]"
        exit 1
    endif

    set pe_input = "$2"
    set pe_outdir = "$3"
    set pe_nfiles = 0
    set pe_beam = 10.6041
    set pe_run = 0
    set pe_qadb = 0
    set pe_ismc = 0
    set pe_workers = 1
    set pe_mx2min = -1.0
    set pe_mx2max = 2.0
    set pe_keep_txt = 0

    if ( $#argv >= 4 ) set pe_nfiles = "$4"
    if ( $#argv >= 5 ) set pe_beam = "$5"
    if ( $#argv >= 6 ) set pe_run = "$6"
    if ( $#argv >= 7 ) set pe_qadb = "$7"
    if ( $#argv >= 8 ) set pe_ismc = "$8"
    if ( $#argv >= 9 ) set pe_workers = "$9"
    if ( $#argv >= 10 ) set pe_mx2min = "$argv[10]"
    if ( $#argv >= 11 ) set pe_mx2max = "$argv[11]"
    if ( $#argv >= 12 ) set pe_keep_txt = "$argv[12]"

    echo "Pulling the latest changes from the repository..."
    git pull
    echo "Sourcing qadb..."
    module load qadb/3.4.1

    ./processing_scripts/run_photon_efficiency_parallel.sh "$pe_input" "$pe_outdir" "$pe_nfiles" "$pe_beam" "$pe_run" "$pe_qadb" "$pe_ismc" "$pe_workers" "$pe_mx2min" "$pe_mx2max" "$pe_keep_txt"
    exit $status
endif

# Set the first argument to be processing two particles if not provided
if ( $#argv < 1 ) then
    set arg1 = "processing_scripts/processing_two_particles.groovy"
    echo "Warning: First argument not provided. Using default: $arg1, which processes two particle events."
else
    set arg1 = "$1"
endif
# Determine the third argument for ./convert_txt_to_root based on arg1
# Initialize to 0 as default value
set convert_arg3 = 0

# Set convert_arg3 based on the value of arg1
if ($arg1 == "processing_scripts/processing_inclusive.groovy") then
    set convert_arg3 = 0
else if ($arg1 == "processing_scripts/processing_mc_inclusive.groovy") then
    set convert_arg3 = 0
else if ($arg1 == "processing_scripts/processing_two_particles.groovy") then
    set convert_arg3 = 1
else if ($arg1 == "processing_scripts/processing_mc_two_particles.groovy") then
    set convert_arg3 = 1
else if ($arg1 == "processing_scripts/processing_three_particles.groovy") then
    set convert_arg3 = 2
else if ($arg1 == "processing_scripts/processing_mc_three_particles.groovy") then
    set convert_arg3 = 2
else if ($arg1 == "processing_scripts/processing_four_particles.groovy") then
    set convert_arg3 = 3
else if ($arg1 == "processing_scripts/processing_dvcs.groovy") then
    set convert_arg3 = 4 # dvcs
else if ($arg1 == "processing_scripts/processing_mc_dvcs.groovy") then
    set convert_arg3 = 4 # dvcs
else if ($arg1 == "processing_scripts/processing_exclusive_pi0.groovy") then
    set convert_arg3 = 5 # eppi0
else if ($arg1 == "processing_scripts/processing_calibration.groovy") then
    set convert_arg3 = 6 # calibration
else if ($arg1 == "processing_scripts/processing_dvcs_calibration.groovy") then
    set convert_arg3 = 6 # calibration
else
    echo "Error: unrecognized processing script: $arg1"
    exit 1
endif

# determine if Monte Carlo
set is_mc = 0;
if ($arg1 == "processing_scripts/processing_mc_inclusive.groovy") then
    set is_mc = 1;
else if ($arg1 == "processing_scripts/processing_mc_two_particles.groovy") then
    set is_mc = 1;
else if ($arg1 == "processing_scripts/processing_mc_three_particles.groovy") then
    set is_mc = 1;
else if ($arg1 == "processing_scripts/processing_mc_dvcs.groovy") then
    set is_mc = 1;
endif

# Set the second argument to default to the RGA Fall2018 inbending nSidis skim if not provided
if ( $#argv < 2 ) then
    set arg2 = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis/"
    echo "Warning: Second argument not provided. Using default: $arg2, the RGAFa18 pass-2 Inbending nSidis skim."
else
    set arg2 = "$2"
endif

echo "Pulling the latest changes from the repository..."
git pull
echo "Sourcing qadb..."
module load qadb/3.4.1

g++ `root-config --cflags --libs` -o processing_scripts/convert_txt_to_root processing_scripts/convert_txt_to_root.cpp

echo "$arg1" "$arg2" "$3.txt" "$4" "$5" "$6" "$7"
# execute command based on number of entries (or dvcs/eppi0/calibration designation)
if ($arg1 == "processing_scripts/processing_inclusive.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3.txt" "$4" "$5" "$6" "$7" 
    # Run the convert_txt_to_root program
    set txt_file = "$3.txt"
    set root_file = "$3.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_mc_inclusive.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3.txt" "$4" "$5" "$6" "$7"
    # Run the convert_txt_to_root program
    set txt_file = "$3.txt"
    set root_file = "$3.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_two_particles.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3" "$4.txt" "$5" "$6" "$7" "$8"
    # Run the convert_txt_to_root program
    set txt_file = "$4.txt"
    set root_file = "$4.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_mc_two_particles.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3" "$4.txt" "$5" "$6" "$7" "$8" "$9"
    # Run the convert_txt_to_root program
    set txt_file = "$4.txt"
    set root_file = "$4.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_three_particles.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3" "$4" "$5.txt" "$6" "$7" "$8" "$9"
    # Run the convert_txt_to_root program
    set txt_file = "$5.txt"
    set root_file = "$5.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_mc_three_particles.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3" "$4" "$5.txt" "$6" "$7" "$8" "$9"
    # Run the convert_txt_to_root program
    set txt_file = "$5.txt"
    set root_file = "$5.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_four_particles.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3" "$4" "$5" "$6.txt" "$7" "$8"
    # Run the convert_txt_to_root program
    set txt_file = "$6.txt"
    set root_file = "$6.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_dvcs.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3.txt" "$4" "$5" "$6" "$7"
    # Run the convert_txt_to_root program
    set txt_file = "$3.txt"
    set root_file = "$3.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_mc_dvcs.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3.txt" "$4" "$5" "$6" "$7"
    # Run the convert_txt_to_root program
    set txt_file = "$3.txt"
    set root_file = "$3.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_exclusive_pi0.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3.txt" "$4" "$5" "$6" "$7"
    # Run the convert_txt_to_root program
    set txt_file = "$3.txt"
    set root_file = "$3.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_calibration.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3.txt" "$4" "$5" "$6" "$7"
    # Run the convert_txt_to_root program
    set txt_file = "$3.txt"
    set root_file = "$3.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
else if ($arg1 == "processing_scripts/processing_dvcs_calibration.groovy") then
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar "$arg1" "$arg2" "$3.txt" "$4" "$5" "$6" "$7"
    # Run the convert_txt_to_root program
    set txt_file = "$3.txt"
    set root_file = "$3.root"
    ./processing_scripts/convert_txt_to_root $txt_file $root_file $convert_arg3 $is_mc
endif
