#include "python_exclusivity_runner.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/wait.h>
#endif

namespace {

namespace fs = std::filesystem;
using json = nlohmann::json;

std::string shell_quote(const std::string& value) {
    std::string quoted = "'";
    for (const char ch : value) {
        if (ch == '\'') {
            quoted += "'\\''";
        } else {
            quoted += ch;
        }
    }
    quoted += "'";
    return quoted;
}

std::string format_decimal(const double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(8) << value;
    std::string result = stream.str();
    while (!result.empty() && result.back() == '0') result.pop_back();
    if (!result.empty() && result.back() == '.') result.pop_back();
    return result;
}

bool is_regular_nonempty_file(const fs::path& path) {
    std::error_code error;
    if (!fs::is_regular_file(path, error) || error) return false;
    const auto size = fs::file_size(path, error);
    return !error && size > 0;
}

bool validate_json_object(const fs::path& path, const std::string& label) {
    if (!is_regular_nonempty_file(path)) {
        std::cerr << "[python_exclusivity_runner] ERROR: Missing or empty "
                  << label << ": " << path << '\n';
        return false;
    }

    try {
        std::ifstream input(path);
        if (!input) {
            std::cerr << "[python_exclusivity_runner] ERROR: Cannot open "
                      << label << ": " << path << '\n';
            return false;
        }

        json document;
        input >> document;
        if (!document.is_object() || document.empty()) {
            std::cerr << "[python_exclusivity_runner] ERROR: " << label
                      << " is not a nonempty JSON object: " << path << '\n';
            return false;
        }
    } catch (const std::exception& error) {
        std::cerr << "[python_exclusivity_runner] ERROR: Invalid JSON in "
                  << label << " (" << path << "): " << error.what() << '\n';
        return false;
    }

    return true;
}

bool atomic_copy(const fs::path& source, const fs::path& destination) {
    std::error_code error;
    fs::create_directories(destination.parent_path(), error);
    if (error) {
        std::cerr << "[python_exclusivity_runner] ERROR: Cannot create directory "
                  << destination.parent_path() << ": " << error.message() << '\n';
        return false;
    }

    fs::path temporary = destination;
    temporary += ".tmp";

    fs::remove(temporary, error);
    error.clear();
    fs::copy_file(source, temporary, fs::copy_options::overwrite_existing, error);
    if (error) {
        std::cerr << "[python_exclusivity_runner] ERROR: Cannot stage "
                  << source << " as " << temporary << ": "
                  << error.message() << '\n';
        return false;
    }

    fs::rename(temporary, destination, error);
    if (error) {
        // Some filesystems do not replace an existing target on rename.
        std::error_code remove_error;
        fs::remove(destination, remove_error);
        error.clear();
        fs::rename(temporary, destination, error);
    }

    if (error) {
        std::cerr << "[python_exclusivity_runner] ERROR: Cannot install "
                  << destination << ": " << error.message() << '\n';
        std::error_code cleanup_error;
        fs::remove(temporary, cleanup_error);
        return false;
    }

    return true;
}

int normalized_exit_code(const int system_status) {
    if (system_status == -1) return -1;
#if defined(__unix__) || defined(__APPLE__)
    if (WIFEXITED(system_status)) return WEXITSTATUS(system_status);
    if (WIFSIGNALED(system_status)) return 128 + WTERMSIG(system_status);
#endif
    return system_status;
}

bool validate_options(const PythonExclusivityOptions& options) {
    if (options.python_executable.empty()) {
        std::cerr << "[python_exclusivity_runner] ERROR: python_executable is empty.\n";
        return false;
    }
    if (options.script_path.empty()) {
        std::cerr << "[python_exclusivity_runner] ERROR: script_path is empty.\n";
        return false;
    }
    if (options.global_cuts_json.empty()) {
        std::cerr << "[python_exclusivity_runner] ERROR: global_cuts_json is empty.\n";
        return false;
    }
    if (options.output_directory.empty() || options.install_directory.empty()) {
        std::cerr << "[python_exclusivity_runner] ERROR: Output directory is empty.\n";
        return false;
    }
    if (options.workers < 1 || options.workers > 7) {
        std::cerr << "[python_exclusivity_runner] ERROR: workers must be in [1, 7], got "
                  << options.workers << ".\n";
        return false;
    }

    const bool finite =
        std::isfinite(options.tight_containment) &&
        std::isfinite(options.nominal_containment) &&
        std::isfinite(options.loose_containment);
    if (!finite ||
        !(0.0 < options.tight_containment &&
          options.tight_containment < options.nominal_containment &&
          options.nominal_containment < options.loose_containment &&
          options.loose_containment < 1.0)) {
        std::cerr
            << "[python_exclusivity_runner] ERROR: Required containment ordering is "
            << "0 < tight < nominal < loose < 1. Received tight="
            << options.tight_containment << ", nominal="
            << options.nominal_containment << ", loose="
            << options.loose_containment << ".\n";
        return false;
    }

    return true;
}

bool all_installed_outputs_are_valid(const fs::path& install_directory) {
    return
        validate_json_object(install_directory / "combined_cuts_90.json", "installed tight cuts") &&
        validate_json_object(install_directory / "combined_cuts_95.json", "installed nominal cuts") &&
        validate_json_object(install_directory / "combined_cuts_98.json", "installed loose cuts") &&
        validate_json_object(install_directory / "combined_cuts.json", "installed compatibility cuts");
}

} // namespace

bool run_python_exclusivity_analysis(
    const PythonExclusivityOptions& options) {

    std::cout << "\n============================================================\n"
              << "Python exclusivity template-fit analysis\n"
              << "============================================================\n";

    if (!options.enabled) {
        std::cout << "[python_exclusivity_runner] Stage disabled; validating existing outputs.\n";
        return all_installed_outputs_are_valid(options.install_directory);
    }

    if (!validate_options(options)) return false;

    const fs::path script_path(options.script_path);
    const fs::path global_cuts_path(options.global_cuts_json);
    const fs::path output_directory(options.output_directory);
    const fs::path install_directory(options.install_directory);

    if (!is_regular_nonempty_file(script_path)) {
        std::cerr << "[python_exclusivity_runner] ERROR: Python script not found: "
                  << script_path << '\n';
        return false;
    }
    if (!is_regular_nonempty_file(global_cuts_path)) {
        std::cerr << "[python_exclusivity_runner] ERROR: Global-cut JSON not found: "
                  << global_cuts_path << '\n';
        return false;
    }

    if (!options.force_rerun && all_installed_outputs_are_valid(install_directory)) {
        std::cout << "[python_exclusivity_runner] Reusing valid installed cut JSONs.\n";
        return true;
    }

    std::error_code error;
    fs::create_directories(output_directory, error);
    if (error) {
        std::cerr << "[python_exclusivity_runner] ERROR: Cannot create "
                  << output_directory << ": " << error.message() << '\n';
        return false;
    }
    error.clear();
    fs::create_directories(install_directory, error);
    if (error) {
        std::cerr << "[python_exclusivity_runner] ERROR: Cannot create "
                  << install_directory << ": " << error.message() << '\n';
        return false;
    }

    std::cout << "Script:              " << script_path << '\n'
              << "Global cuts:         " << global_cuts_path << '\n'
              << "Fit output:          " << output_directory << '\n'
              << "Canonical JSON dir:  " << install_directory << '\n'
              << "Tight containment:   " << 100.0 * options.tight_containment << "%\n"
              << "Nominal containment: " << 100.0 * options.nominal_containment << "%\n"
              << "Loose containment:   " << 100.0 * options.loose_containment << "%\n"
              << "Workers:             " << options.workers << '\n';

    std::ostringstream command;
    command
        << shell_quote(options.python_executable) << ' '
        << shell_quote(script_path.string())
        << " --global-cuts-json " << shell_quote(global_cuts_path.string())
        << " --output-dir " << shell_quote(output_directory.string())
        << " --tight-containment " << format_decimal(options.tight_containment)
        << " --nominal-containment " << format_decimal(options.nominal_containment)
        << " --loose-containment " << format_decimal(options.loose_containment)
        << " --workers " << options.workers;

    std::cout << "[python_exclusivity_runner] Executing:\n"
              << command.str() << "\n\n" << std::flush;

    const int raw_status = std::system(command.str().c_str());
    const int exit_code = normalized_exit_code(raw_status);
    if (exit_code != 0) {
        std::cerr << "[python_exclusivity_runner] ERROR: Python stage exited with code "
                  << exit_code << ".\n";
        return false;
    }

    // This path matches the production Python script's current output hierarchy.
    const fs::path generated_json_directory =
        output_directory / "iterative_cuts" / "jsons";
    const fs::path generated_tight = generated_json_directory / "combined_cuts_90.json";
    const fs::path generated_nominal = generated_json_directory / "combined_cuts_95.json";
    const fs::path generated_loose = generated_json_directory / "combined_cuts_98.json";

    if (!validate_json_object(generated_tight, "generated 90% cuts") ||
        !validate_json_object(generated_nominal, "generated 95% cuts") ||
        !validate_json_object(generated_loose, "generated 98% cuts")) {
        return false;
    }

    if (!atomic_copy(generated_tight, install_directory / "combined_cuts_90.json") ||
        !atomic_copy(generated_nominal, install_directory / "combined_cuts_95.json") ||
        !atomic_copy(generated_loose, install_directory / "combined_cuts_98.json") ||
        !atomic_copy(generated_nominal, install_directory / "combined_cuts.json")) {
        return false;
    }

    if (!all_installed_outputs_are_valid(install_directory)) return false;

    std::cout << "[python_exclusivity_runner] Installed and validated:\n"
              << "  " << install_directory / "combined_cuts_90.json" << '\n'
              << "  " << install_directory / "combined_cuts_95.json" << '\n'
              << "  " << install_directory / "combined_cuts_98.json" << '\n'
              << "  " << install_directory / "combined_cuts.json" << '\n'
              << "============================================================\n";

    return true;
}
