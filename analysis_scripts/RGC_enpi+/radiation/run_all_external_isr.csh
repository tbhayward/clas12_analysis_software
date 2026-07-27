#!/bin/tcsh -f

# =============================================================================
# Apply the external-ISR transformation to every momentum-corrected RGC ISR
# ROOT file produced by run_all_momentum_corrections.csh.
#
# Input selection:
#   *_ISR_mom_corrections.root
#
# Excluded:
#   * Files that already contain "externalISR" in the basename
#
# Output naming is handled by apply_external_isr.py:
#   *_ISR_mom_corrections.root
#       -> *_ISR_externalISR_mom_corrections.root
#
# The transformed ROOT files, provenance JSON files, and validation directories
# are written into the same paper_versions directory as the inputs.
# =============================================================================

# Keep this launcher, apply_external_isr.py, and the geometry JSON together.
set script_dir = "${0:h}"
if ("${script_dir}" == "") then
    set script_dir = "."
endif
set script_dir = "`cd "${script_dir}" && pwd`"

set data_dir = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions"

set external_isr_script = "${script_dir}/apply_external_isr.py"
set geometry_json = "${script_dir}/external_radiation_geometry.json"

if (! -f "${external_isr_script}") then
    echo "FATAL: External-ISR script does not exist:"
    echo "       ${external_isr_script}"
    exit 1
endif

if (! -f "${geometry_json}") then
    echo "FATAL: Geometry JSON does not exist:"
    echo "       ${geometry_json}"
    exit 1
endif

if (! -d "${data_dir}") then
    echo "FATAL: Data directory does not exist:"
    echo "       ${data_dir}"
    exit 1
endif

set log_dir = "${data_dir}/external_isr_logs"
mkdir -p "${log_dir}"

if ($status != 0) then
    echo "FATAL: Could not create log directory:"
    echo "       ${log_dir}"
    exit 1
endif

set n_found = 0
set n_succeeded = 0
set n_failed = 0
set failed_files = ()

echo "======================================================================"
echo "RGC external-ISR batch"
echo "======================================================================"
echo "Input/output directory: ${data_dir}"
echo "Python script:          ${external_isr_script}"
echo "Geometry JSON:          ${geometry_json}"
echo "Log directory:          ${log_dir}"
echo "======================================================================"
echo ""

foreach input_file ("${data_dir}"/*_ISR_mom_corrections.root)

    # tcsh can preserve an unmatched glob literally, so guard against that.
    if (! -f "${input_file}") then
        continue
    endif

    set input_name = "${input_file:t}"

    # Do not recursively process files already transformed by this stage.
    if ("${input_name}" =~ "*externalISR*") then
        echo "Skipping previously transformed file: ${input_name}"
        continue
    endif

    @ n_found++

    set input_stem = "${input_name:r}"
    set output_stem = "`echo "${input_stem}" | sed 's/_ISR_/_ISR_externalISR_/'`"
    set output_file = "${data_dir}/${output_stem}.root"
    set provenance_file = "${data_dir}/${output_stem}.provenance.json"
    set validation_dir = "${data_dir}/${output_stem}_validation"
    set log_file = "${log_dir}/${output_stem}.log"

    echo ""
    echo "----------------------------------------------------------------------"
    echo "[${n_found}] Processing"
    echo "Input:      ${input_file}"
    echo "Output:     ${output_file}"
    echo "Provenance: ${provenance_file}"
    echo "Validation: ${validation_dir}"
    echo "Log:        ${log_file}"
    echo "----------------------------------------------------------------------"

    # apply_external_isr.py intentionally refuses to overwrite existing ROOT
    # or provenance files. Remove all products from a previous attempt so a
    # rerun is deterministic and complete.
    if (-e "${output_file}") then
        rm -f "${output_file}"
    endif
    if (-e "${provenance_file}") then
        rm -f "${provenance_file}"
    endif
    if (-d "${validation_dir}") then
        rm -rf "${validation_dir}"
    endif

    python3 "${external_isr_script}" \
        "${input_file}" \
        --geometry-json "${geometry_json}" \
        --output-file "${output_file}" \
        --validation-dir "${validation_dir}" \
        >&! "${log_file}"

    set run_status = $status

    if (${run_status} == 0) then
        @ n_succeeded++
        echo "SUCCESS: ${output_file:t}"
    else
        @ n_failed++
        set failed_files = (${failed_files} "${input_name}")
        echo "FAILED: ${input_name}"
        echo "See log: ${log_file}"

        # Remove incomplete products from a failed transformation.
        if (-e "${output_file}") then
            rm -f "${output_file}"
        endif
        if (-e "${provenance_file}") then
            rm -f "${provenance_file}"
        endif
        if (-d "${validation_dir}") then
            rm -rf "${validation_dir}"
        endif
    endif

end

echo ""
echo "======================================================================"
echo "Batch complete"
echo "======================================================================"
echo "Eligible ISR files found: ${n_found}"
echo "Succeeded:                ${n_succeeded}"
echo "Failed:                   ${n_failed}"
echo "Output directory:         ${data_dir}"
echo "Log directory:            ${log_dir}"

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
    echo "WARNING: No *_ISR_mom_corrections.root files were found."
    exit 2
endif

exit 0
