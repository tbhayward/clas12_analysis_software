#include "external_scripts_runner.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/wait.h>
#endif

namespace fs = std::filesystem;

namespace {

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

bool command_succeeded(const int status, int& exit_code) {
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

bool run_one_script(
    const fs::path& scripts_directory,
    const fs::path& analysis_csv,
    const std::string& python_executable,
    const std::string& script_name,
    const std::string& extra_arguments
) {
    const fs::path script_path = scripts_directory / script_name;

    if (!fs::is_regular_file(script_path)) {
        std::cerr << "[external scripts] ERROR: missing script: "
                  << script_path << '\n';
        return false;
    }

    std::string command =
        "cd " + shell_quote(scripts_directory.string()) +
        " && " + shell_quote(python_executable) +
        " " + shell_quote(script_name) +
        " " + shell_quote(analysis_csv.string());

    if (!extra_arguments.empty()) {
        command += " " + extra_arguments;
    }

    std::cout << "[external scripts] Running " << script_name << "...\n";
    std::cout.flush();

    const int status = std::system(command.c_str());
    int exit_code = -1;

    if (!command_succeeded(status, exit_code)) {
        std::cerr << "[external scripts] ERROR: " << script_name
                  << " failed with exit code " << exit_code << ".\n";
        return false;
    }

    std::cout << "[external scripts] Finished " << script_name << ".\n";
    return true;
}

} // namespace

bool run_external_cross_section_scripts(
    const std::string& analysis_csv,
    const ExternalScriptOptions& options
) {
    std::error_code ec;

    const fs::path csv_absolute = fs::absolute(analysis_csv, ec);
    if (ec || !fs::is_regular_file(csv_absolute)) {
        std::cerr << "[external scripts] ERROR: analysis CSV does not exist: "
                  << analysis_csv << '\n';
        return false;
    }

    ec.clear();
    const fs::path scripts_absolute =
        fs::absolute(options.scripts_directory, ec);
    if (ec || !fs::is_directory(scripts_absolute)) {
        std::cerr << "[external scripts] ERROR: scripts directory does not exist: "
                  << options.scripts_directory << '\n';
        return false;
    }

    const std::string integrated_arguments =
        options.include_bin_to_bin_systematics
            ? "--include-bin-to-bin-sys"
            : "";

    const std::string clas6_arguments =
        options.use_simple_clas6_cross_check
            ? "--simple"
            : "";

    if (!run_one_script(
            scripts_absolute,
            csv_absolute,
            options.python_executable,
            "integrated_cross_section_study.py",
            integrated_arguments)) {
        return false;
    }

    if (!run_one_script(
            scripts_absolute,
            csv_absolute,
            options.python_executable,
            "clas6_cross_check.py",
            clas6_arguments)) {
        return false;
    }

    if (!run_one_script(
            scripts_absolute,
            csv_absolute,
            options.python_executable,
            "hall_a_cross_check.py",
            "")) {
        return false;
    }

    return true;
}
