#include "systematics_runner.h"

#include <cstdlib>
#include <chrono>
#include <iomanip>
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

} // namespace

bool run_main_systematics(
    const std::string& analysis_csv,
    const SystematicsRunnerOptions& options
) {
    std::error_code ec;

    const fs::path executable_absolute = fs::absolute(options.executable, ec);
    if (ec || !fs::is_regular_file(executable_absolute)) {
        std::cerr << "[systematics runner] ERROR: systematics executable does not exist: "
                  << options.executable << '\n';
        return false;
    }

    ec.clear();
    const fs::path csv_absolute = fs::absolute(analysis_csv, ec);
    if (ec || !fs::is_regular_file(csv_absolute)) {
        std::cerr << "[systematics runner] ERROR: analysis CSV does not exist: "
                  << analysis_csv << '\n';
        return false;
    }

    ec.clear();
    const fs::path pass1_absolute =
        fs::absolute(options.pass1_systematics_csv, ec);
    if (ec || !fs::is_regular_file(pass1_absolute)) {
        std::cerr << "[systematics runner] ERROR: pass-1 systematics CSV does not exist: "
                  << options.pass1_systematics_csv << '\n';
        return false;
    }

    const std::string command =
        shell_quote(executable_absolute.string()) + " " +
        shell_quote(csv_absolute.string()) + " " +
        shell_quote(pass1_absolute.string());

    std::cout << "[systematics runner] Running main_systematics on the freshly generated CSV...\n";
    std::cout.flush();
    const auto t0 = std::chrono::steady_clock::now();

    const int status = std::system(command.c_str());
    int exit_code = -1;

    if (!command_succeeded(status, exit_code)) {
        std::cerr << "[systematics runner] ERROR: main_systematics failed with exit code "
                  << exit_code << ".\n";
        return false;
    }

    const double elapsed_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    std::cout << "[systematics runner] main_systematics finished successfully in "
              << std::fixed << std::setprecision(1) << elapsed_s
              << " s." << std::defaultfloat << std::setprecision(6) << "\n";
    return true;
}
