#!/bin/tcsh -f

# =============================================================================
# Apply the finalized RGC electron and pi+ momentum corrections to every
# eligible RGC ROOT file.
#
# Excluded:
#   * Files ending in _2.root
#   * temp_mc.root
#   * RGA files, because apply_momentum_corrections.py only supports RGC runs
#
# Output naming:
#   input_name.root -> input_name_mom_corrections.root
# =============================================================================

set script_dir = "/u/home/thayward/clas12_analysis_software/analysis_scripts/RGC_enpi+/momentum_corrections"

set input_dir = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+"

set output_dir = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions"

set correction_script = "${script_dir}/apply_momentum_corrections.py"

# Move to the correction-script directory so its default JSON paths resolve
# correctly:
#
#   output/electron_diagnostics/json/electron_correction_models.json
#   output/piplus_diagnostics/json/pion_correction_models.json

cd "${script_dir}"

if ($status != 0) then
    echo "FATAL: Could not enter ${script_dir}"
    exit 1
endif

if (! -f "${correction_script}") then
    echo "FATAL: Correction script does not exist:"
    echo "       ${correction_script}"
    exit 1
endif

if (! -d "${input_dir}") then
    echo "FATAL: Input directory does not exist:"
    echo "       ${input_dir}"
    exit 1
endif

mkdir -p "${output_dir}"

if ($status != 0) then
    echo "FATAL: Could not create output directory:"
    echo "       ${output_dir}"
    exit 1
endif

set log_dir = "${output_dir}/momentum_correction_logs"
mkdir -p "${log_dir}"

set n_found = 0
set n_succeeded = 0
set n_failed = 0
set failed_files = ()

echo "======================================================================"
echo "RGC momentum-correction batch"
echo "======================================================================"
echo "Input directory:  ${input_dir}"
echo "Output directory: ${output_dir}"
echo "Python script:    ${correction_script}"
echo "======================================================================"
echo ""

foreach input_file ("${input_dir}"/rgc_*.root)

    # tcsh can preserve an unmatched glob literally, so guard against that.
    if (! -f "${input_file}") then
        continue
    endif

    set input_name = "${input_file:t}"

    # Exclude files whose basename ends in _2.root.
    if ("${input_name}" =~ "*_2.root") then
        echo "Skipping secondary file: ${input_name}"
        continue
    endif

    # Exclude temp_mc.root explicitly. This is redundant with the rgc_*.root
    # glob, but retained here to make the intended selection unambiguous.
    if ("${input_name}" == "temp_mc.root") then
        echo "Skipping temporary file: ${input_name}"
        continue
    endif

    @ n_found++

    set name_without_extension = "${input_name:r}"
    set output_name = "${name_without_extension}_mom_corrections.root"
    set output_file = "${output_dir}/${output_name}"
    set log_file = "${log_dir}/${name_without_extension}_mom_corrections.log"

    echo ""
    echo "----------------------------------------------------------------------"
    echo "[${n_found}] Processing"
    echo "Input:  ${input_file}"
    echo "Output: ${output_file}"
    echo "Log:    ${log_file}"
    echo "----------------------------------------------------------------------"

    python3 "${correction_script}" \
        "${input_file}" \
        "${output_file}" \
        --overwrite \
        >&! "${log_file}"

    set run_status = $status

    if (${run_status} == 0) then
        @ n_succeeded++
        echo "SUCCESS: ${output_name}"
    else
        @ n_failed++
        set failed_files = (${failed_files} "${input_name}")
        echo "FAILED: ${input_name}"
        echo "See log: ${log_file}"

        # Remove an incomplete output file if the Python process failed.
        if (-e "${output_file}") then
            rm -f "${output_file}"
        endif
    endif

end

echo ""
echo "======================================================================"
echo "Batch complete"
echo "======================================================================"
echo "Eligible files found: ${n_found}"
echo "Succeeded:            ${n_succeeded}"
echo "Failed:               ${n_failed}"
echo "Output directory:     ${output_dir}"
echo "Log directory:        ${log_dir}"

if (${n_failed} > 0) then
    echo ""
    echo "Failed input files:"

    foreach failed_file (${failed_files})
        echo "  ${failed_file}"
    end

    exit 1
endif

if (${n_found} == 0) then
    echo ""
    echo "WARNING: No eligible RGC ROOT files were found."
    exit 2
endif

exit 0