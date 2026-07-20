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

namespace fs = std::filesystem;

namespace {

constexpr const char* kLogPrefix = "[python exclusivity]";

std::string shell_quote(const std::string& text) {
    std::string quoted;
    quoted.reserve(text.size() + 2);
    quoted.push_back('\'');

    for (const char ch : text) {
        if (ch == '\'') {
            quoted += "'\\''";
        } else {
            quoted.push_back(ch);
        }
    }

    quoted.push_back('\'');
    return quoted;
}

bool decode_command_status(const int status, int& exit_code) {
    exit_code = -1;

    if (status == -1) {
        return false;
    }

#if defined(__unix__) || defined(__APPLE__)
    if (WIFEXITED(status)) {
        exit_code = WEXITSTATUS(status);
        return exit_code == 0;
    }

    if (WIFSIGNALED(status)) {
        exit_code = 128 + WTERMSIG(status);
    }

    return false;
#else
    exit_code = status;
    return status == 0;
#endif
}

int execute_command(const std::string& command) {
    std::cout << kLogPrefix << " Command:\n" << command << "\n";
    std::cout.flush();

    const int status = std::system(command.c_str());
    int exit_code = -1;
    decode_command_status(status, exit_code);
    return exit_code;
}

bool validate_containments(const PythonExclusivityOptions& options) {
    const auto valid = [](const double value) {
        return std::isfinite(value) && value > 0.0 && value < 1.0;
    };

    if (!valid(options.tight_containment) ||
        !valid(options.nominal_containment) ||
        !valid(options.loose_containment)) {
        std::cerr << kLogPrefix
                  << " ERROR: all containment values must be finite and strictly "
                     "between zero and one.\n";
        return false;
    }

    if (!(options.tight_containment < options.nominal_containment &&
          options.nominal_containment < options.loose_containment)) {
        std::cerr << kLogPrefix
                  << " ERROR: containment values must satisfy tight < nominal < loose.\n";
        return false;
    }

    if (options.workers < 1) {
        std::cerr << kLogPrefix << " ERROR: workers must be at least one.\n";
        return false;
    }

    return true;
}

bool validate_json_file(const fs::path& path) {
    std::error_code ec;
    if (!fs::is_regular_file(path, ec) || ec) {
        std::cerr << kLogPrefix << " ERROR: missing JSON file: " << path << '\n';
        return false;
    }

    const std::uintmax_t size = fs::file_size(path, ec);
    if (ec || size == 0) {
        std::cerr << kLogPrefix << " ERROR: JSON file is empty or unreadable: "
                  << path << '\n';
        return false;
    }

    try {
        std::ifstream input(path);
        if (!input) {
            std::cerr << kLogPrefix << " ERROR: cannot open JSON file: "
                      << path << '\n';
            return false;
        }

        nlohmann::json document;
        input >> document;

        if (!document.is_object() || document.empty()) {
            std::cerr << kLogPrefix
                      << " ERROR: expected a nonempty top-level JSON object: "
                      << path << '\n';
            return false;
        }
    } catch (const std::exception& error) {
        std::cerr << kLogPrefix << " ERROR: malformed JSON file " << path
                  << ": " << error.what() << '\n';
        return false;
    }

    return true;
}

bool copy_file_atomically(const fs::path& source, const fs::path& destination) {
    std::error_code ec;
    fs::create_directories(destination.parent_path(), ec);
    if (ec) {
        std::cerr << kLogPrefix << " ERROR: cannot create directory "
                  << destination.parent_path() << ": " << ec.message() << '\n';
        return false;
    }

    fs::path temporary = destination;
    temporary += ".tmp";

    fs::remove(temporary, ec);
    ec.clear();

    fs::copy_file(source, temporary, fs::copy_options::overwrite_existing, ec);
    if (ec) {
        std::cerr << kLogPrefix << " ERROR: cannot copy " << source << " to "
                  << temporary << ": " << ec.message() << '\n';
        return false;
    }

    fs::rename(temporary, destination, ec);
    if (ec) {
        // Some filesystems do not replace an existing destination on rename.
        std::error_code remove_ec;
        fs::remove(destination, remove_ec);
        ec.clear();
        fs::rename(temporary, destination, ec);
    }

    if (ec) {
        std::cerr << kLogPrefix << " ERROR: cannot install " << destination
                  << ": " << ec.message() << '\n';
        std::error_code cleanup_ec;
        fs::remove(temporary, cleanup_ec);
        return false;
    }

    return true;
}

bool installed_outputs_are_valid(const fs::path& install_directory) {
    return validate_json_file(install_directory / "combined_cuts_90.json") &&
           validate_json_file(install_directory / "combined_cuts_95.json") &&
           validate_json_file(install_directory / "combined_cuts_98.json") &&
           validate_json_file(install_directory / "combined_cuts.json");
}

std::string format_double(const double value) {
    std::ostringstream output;
    output << std::fixed << std::setprecision(6) << value;
    return output.str();
}

} // namespace

bool run_python_exclusivity_analysis(
    const PythonExclusivityOptions& options
) {
    if (!options.enabled) {
        std::cout << kLogPrefix << " Stage disabled; using existing cut JSONs.\n";
        return installed_outputs_are_valid(options.install_directory);
    }

    if (!validate_containments(options)) {
        return false;
    }

    std::error_code ec;

    const fs::path script_absolute = fs::absolute(options.script_path, ec);
    if (ec || !fs::is_regular_file(script_absolute)) {
        std::cerr << kLogPrefix << " ERROR: Python script does not exist: "
                  << options.script_path << '\n';
        return false;
    }

    ec.clear();
    const fs::path global_cuts_absolute =
        fs::absolute(options.global_cuts_json, ec);
    if (ec || !fs::is_regular_file(global_cuts_absolute)) {
        std::cerr << kLogPrefix << " ERROR: global-cuts JSON does not exist: "
                  << options.global_cuts_json << '\n';
        return false;
    }

    ec.clear();
    const fs::path output_absolute =
        fs::absolute(options.output_directory, ec);
    if (ec) {
        std::cerr << kLogPrefix << " ERROR: cannot resolve output directory: "
                  << options.output_directory << '\n';
        return false;
    }

    ec.clear();
    const fs::path install_absolute =
        fs::absolute(options.install_directory, ec);
    if (ec) {
        std::cerr << kLogPrefix << " ERROR: cannot resolve install directory: "
                  << options.install_directory << '\n';
        return false;
    }

    if (!options.force_rerun && installed_outputs_are_valid(install_absolute)) {
        std::cout << kLogPrefix
                  << " Existing validated cut JSONs found; skipping Python rerun.\n";
        return true;
    }

    fs::create_directories(output_absolute, ec);
    if (ec) {
        std::cerr << kLogPrefix << " ERROR: cannot create output directory "
                  << output_absolute << ": " << ec.message() << '\n';
        return false;
    }

    ec.clear();
    fs::create_directories(install_absolute, ec);
    if (ec) {
        std::cerr << kLogPrefix << " ERROR: cannot create install directory "
                  << install_absolute << ": " << ec.message() << '\n';
        return false;
    }

    std::cout << "\n============================================================\n"
              << "Python exclusivity optimization\n"
              << "============================================================\n"
              << "Script:        " << script_absolute << '\n'
              << "Global cuts:   " << global_cuts_absolute << '\n'
              << "Output:        " << output_absolute << '\n'
              << "Install:       " << install_absolute << '\n'
              << "Workers:       " << options.workers << '\n'
              << "Tight:         " << 100.0 * options.tight_containment << "%\n"
              << "Nominal:       " << 100.0 * options.nominal_containment << "%\n"
              << "Loose:         " << 100.0 * options.loose_containment << "%\n"
              << "============================================================\n";

    std::string command =
        shell_quote(options.python_executable) +
        " " + shell_quote(script_absolute.string()) +
        " --global-cuts-json " + shell_quote(global_cuts_absolute.string()) +
        " --output-directory " + shell_quote(output_absolute.string()) +
        " --tight-containment " + format_double(options.tight_containment) +
        " --nominal-containment " + format_double(options.nominal_containment) +
        " --loose-containment " + format_double(options.loose_containment) +
        " --workers " + std::to_string(options.workers);

    const int exit_code = execute_command(command);
    if (exit_code != 0) {
        std::cerr << kLogPrefix
                  << " ERROR: Python exclusivity stage failed with exit code "
                  << exit_code << ".\n";
        return false;
    }

    const fs::path generated_json_directory = output_absolute / "jsons";
    const fs::path generated_90 =
        generated_json_directory / "combined_cuts_90.json";
    const fs::path generated_95 =
        generated_json_directory / "combined_cuts_95.json";
    const fs::path generated_98 =
        generated_json_directory / "combined_cuts_98.json";

    if (!validate_json_file(generated_90) ||
        !validate_json_file(generated_95) ||
        !validate_json_file(generated_98)) {
        std::cerr << kLogPrefix
                  << " ERROR: Python completed but did not produce all required "
                     "validated cut JSONs.\n";
        return false;
    }

    if (!copy_file_atomically(
            generated_90, install_absolute / "combined_cuts_90.json") ||
        !copy_file_atomically(
            generated_95, install_absolute / "combined_cuts_95.json") ||
        !copy_file_atomically(
            generated_98, install_absolute / "combined_cuts_98.json") ||
        !copy_file_atomically(
            generated_95, install_absolute / "combined_cuts.json")) {
        return false;
    }

    if (!installed_outputs_are_valid(install_absolute)) {
        std::cerr << kLogPrefix
                  << " ERROR: one or more installed cut JSONs failed final validation.\n";
        return false;
    }

    std::cout << kLogPrefix << " Verified and installed:\n"
              << "  " << install_absolute / "combined_cuts_90.json" << '\n'
              << "  " << install_absolute / "combined_cuts_95.json" << '\n'
              << "  " << install_absolute / "combined_cuts_98.json" << '\n'
              << "  " << install_absolute / "combined_cuts.json"
              << "  (nominal 95%)\n";

    return true;
}
