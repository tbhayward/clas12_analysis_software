#!/bin/csh

source source_file.txt

# Set the first argument to be processing two particles if not provided
if ( $#argv < 1 ) then
    set arg1 = "processing_scripts/processing_two_particles.groovy"
    echo "Warning: First argument not provided. Using default: $arg1, which processes two particle events."
else
    set arg1 = "$1"
endif
# Determine the third argument for ./convert_txt_to_root based on arg1
set convert_arg3 = 0

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
    set convert_arg3 = 4
else if ($arg1 == "processing_scripts/processing_mc_dvcs.groovy") then
    set convert_arg3 = 4
else if ($arg1 == "processing_scripts/processing_exclusive_pi0.groovy") then
    set convert_arg3 = 5
else if ($arg1 == "processing_scripts/processing_calibration.groovy") then
    set convert_arg3 = 6
else if ($arg1 == "processing_scripts/processing_dvcs_calibration.groovy") then
    set convert_arg3 = 6
else if ($arg1 == "processing_scripts/processing_epgamma.groovy") then
    set convert_arg3 = 9
else if ($arg1 == "processing_scripts/processing_mc_epgamma.groovy") then
    set convert_arg3 = 9
else if ($arg1 == "processing_scripts/processing_epgammagamma.groovy") then
    set convert_arg3 = 10
else if ($arg1 == "processing_scripts/processing_mc_epgammagamma.groovy") then
    set convert_arg3 = 10
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
else if ($arg1 == "processing_scripts/processing_mc_epgamma.groovy") then
    set is_mc = 1;
else if ($arg1 == "processing_scripts/processing_mc_epgammagamma.groovy") then
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

echo "Processing script: $arg1"
echo "Input path: $arg2"
if ( $#argv >= 3 ) then
    echo "Output base: $3"
endif
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
else if ($arg1 == "processing_scripts/processing_epgamma.groovy" || $arg1 == "processing_scripts/processing_mc_epgamma.groovy") then
    # ------------------------------------------------------------------
    # Unified data/MC e'p'gammaX processor.
    #
    # Calling convention:
    #   processing.csh <script> <hipo_dir> <output_base>
    #                  [n_files] [beam_energy] [run_override] [qadb_override]
    #
    # Defaults:
    #   n_files       = 0       (all files)
    #   beam_energy   = 10.604  GeV
    #   run_override  = 0       (use RUN::config)
    #   qadb_override = 1
    # ------------------------------------------------------------------
    if ( $#argv < 3 ) then
        echo "ERROR: processing_epgamma requires at least:"
        echo "  processing.csh <script> <hipo_dir> <output_base>"
        exit 1
    endif

    set photon_n_files = 0
    if ( $#argv >= 4 ) set photon_n_files = "$4"

    set photon_beam_energy = 10.604
    if ( $#argv >= 5 ) set photon_beam_energy = "$5"

    set photon_run_override = 0
    if ( $#argv >= 6 ) set photon_run_override = "$6"

    set photon_qadb_override = 1
    if ( $#argv >= 7 ) set photon_qadb_override = "$7"

    set output_base = "$3"
    set quick_txt = "${output_base}_1M.txt"
    set quick_root = "${output_base}_1M.root"
    set full_txt = "${output_base}.txt"
    set full_root = "${output_base}.root"

    echo "Photon-processing resolved arguments:"
    echo "  n_files       = $photon_n_files (0 = all)"
    echo "  beam_energy   = $photon_beam_energy GeV"
    echo "  run_override  = $photon_run_override"
    echo "  qadb_override = $photon_qadb_override"

    echo "=== epgamma quick pass: stopping after 1,000,000 candidate rows ==="
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar \
        "$arg1" "$arg2" "$quick_txt" \
        "$photon_n_files" "$photon_beam_energy" \
        "$photon_run_override" "$photon_qadb_override" \
        "1000000"
    if ($status != 0) exit $status

    ./processing_scripts/convert_txt_to_root \
        "$quick_txt" "$quick_root" $convert_arg3 $is_mc
    if ($status != 0) exit $status

    echo "=== epgamma full pass: restarting from the beginning for all statistics ==="
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar \
        "$arg1" "$arg2" "$full_txt" \
        "$photon_n_files" "$photon_beam_energy" \
        "$photon_run_override" "$photon_qadb_override" \
        "0"
    if ($status != 0) exit $status

    ./processing_scripts/convert_txt_to_root \
        "$full_txt" "$full_root" $convert_arg3 $is_mc

else if ($arg1 == "processing_scripts/processing_epgammagamma.groovy" || $arg1 == "processing_scripts/processing_mc_epgammagamma.groovy") then
    # Same optional argument ordering/defaults as processing_epgamma.groovy:
    # [n_files] [beam_energy] [run_override] [qadb_override]
    if ( $#argv < 3 ) then
        echo "ERROR: processing_epgammagamma requires at least:"
        echo "  processing.csh <script> <hipo_dir> <output_base>"
        exit 1
    endif

    set photon_n_files = 0
    if ( $#argv >= 4 ) set photon_n_files = "$4"

    set photon_beam_energy = 10.604
    if ( $#argv >= 5 ) set photon_beam_energy = "$5"

    set photon_run_override = 0
    if ( $#argv >= 6 ) set photon_run_override = "$6"

    set photon_qadb_override = 1
    if ( $#argv >= 7 ) set photon_qadb_override = "$7"

    set output_base = "$3"
    set quick_txt = "${output_base}_1M.txt"
    set quick_root = "${output_base}_1M.root"
    set full_txt = "${output_base}.txt"
    set full_root = "${output_base}.root"

    echo "Photon-processing resolved arguments:"
    echo "  n_files       = $photon_n_files (0 = all)"
    echo "  beam_energy   = $photon_beam_energy GeV"
    echo "  run_override  = $photon_run_override"
    echo "  qadb_override = $photon_qadb_override"

    echo "=== epgammagamma quick pass: stopping after 1,000,000 candidate rows ==="
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar \
        "$arg1" "$arg2" "$quick_txt" \
        "$photon_n_files" "$photon_beam_energy" \
        "$photon_run_override" "$photon_qadb_override" \
        "1000000"
    if ($status != 0) exit $status

    ./processing_scripts/convert_txt_to_root \
        "$quick_txt" "$quick_root" $convert_arg3 $is_mc
    if ($status != 0) exit $status

    echo "=== epgammagamma full pass: restarting from the beginning for all statistics ==="
    coatjava/bin/run-groovy -cp processing_classes/dist/processing_classes.jar \
        "$arg1" "$arg2" "$full_txt" \
        "$photon_n_files" "$photon_beam_energy" \
        "$photon_run_override" "$photon_qadb_override" \
        "0"
    if ($status != 0) exit $status

    ./processing_scripts/convert_txt_to_root \
        "$full_txt" "$full_root" $convert_arg3 $is_mc

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
