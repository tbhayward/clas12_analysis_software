// model_predictions.cpp

#include "model_predictions.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <array>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <sys/stat.h>

#ifdef _WIN32
#  define POPEN  _popen
#  define PCLOSE _pclose
#else
#  define POPEN  popen
#  define PCLOSE pclose
#endif

// --------- Defaults you can edit ----------
static const char* kDefaultDvcsgenPath = "/u/home/thayward/dvcsgens/dvcsgen_print";
static const char* kDefaultKm15Cli     = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/km15_cli.py";
// -----------------------------------------

static std::string trim_copy(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return std::string();
    return s.substr(a, b - a + 1);
}

static std::string shell_quote_posix(const std::string& s) {
    // Single-quote for /bin/sh. Safe for paths with spaces or special chars.
    // 'abc' -> 'abc'
    // abc'def -> 'abc'"'"'def'
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('\'');
    for (char c : s) {
        if (c == '\'') {
            out.append("'\"'\"'");
        } else {
            out.push_back(c);
        }
    }
    out.push_back('\'');
    return out;
}

// Find an executable path to dvcsgen directory.
// Priority: paths.dvcsgen_path -> env DVCSGEN_PATH -> default.
static std::string resolve_dvcsgen_dir(const ModelPaths& paths) {
    if (!paths.dvcsgen_path.empty()) return paths.dvcsgen_path;
    const char* envp = std::getenv("DVCSGEN_PATH");
    if (envp && std::strlen(envp) > 0) return std::string(envp);
    return std::string(kDefaultDvcsgenPath);
}

// Find an absolute path to km15_cli.py.
// Priority: paths.km15_cli -> env KM15_CLI -> default.
static std::string resolve_km15_cli(const ModelPaths& paths) {
    if (!paths.km15_cli.empty()) return paths.km15_cli;
    const char* envp = std::getenv("KM15_CLI");
    if (envp && std::strlen(envp) > 0) return std::string(envp);
    return std::string(kDefaultKm15Cli);
}

// Run a shell command, capture stdout. Throws on failure if require_ok is true.
static std::string run_and_capture_stdout(const std::string& cmd, bool require_ok) {
    std::array<char, 4096> buf{};
    std::string out;

    FILE* pipe = POPEN(cmd.c_str(), "r");
    if (!pipe) {
        if (require_ok) throw std::runtime_error("popen failed: " + cmd);
        return std::string();
    }

    while (true) {
        size_t n = std::fread(buf.data(), 1, buf.size(), pipe);
        if (n == 0) break;
        out.append(buf.data(), n);
    }

    int rc = PCLOSE(pipe);
    if (require_ok && rc != 0) {
        std::ostringstream oss;
        oss << "Command failed (rc=" << rc << "): " << cmd << "\nOutput:\n" << out;
        throw std::runtime_error(oss.str());
    }

    return out;
}


static bool executable_exists(const std::string& path) {
    if (path.empty()) return false;
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    return (st.st_mode & S_IXUSR) || (st.st_mode & S_IXGRP) || (st.st_mode & S_IXOTH);
}

static std::string shell_command_stdout_quiet(const std::string& cmd) {
    return run_and_capture_stdout(cmd + " 2>/dev/null", /*require_ok=*/false);
}

static std::string find_python3_on_path() {
#ifndef _WIN32
    const std::string out = shell_command_stdout_quiet("command -v python3");
    const std::string py = trim_copy(out);
    if (!py.empty() && executable_exists(py)) return py;
#endif
    return std::string("python3");
}

static std::string resolve_km15_python() {
#ifndef _WIN32
    const char* py_km15_env = std::getenv("PY_KM15");
    static bool warned_bad_py = false;
    if (py_km15_env && std::strlen(py_km15_env) > 0) {
        const std::string requested(py_km15_env);
        if (executable_exists(requested)) return requested;
        if (!warned_bad_py) {
            std::cerr << "[km15_xs] WARNING: PY_KM15 points to a missing/non-executable file: "
                      << requested << "\n";
            std::cerr << "[km15_xs] Falling back to python3 from PATH, matching the external plotting script behavior.\n";
            warned_bad_py = true;
        }
    }
    return find_python3_on_path();
#else
    return std::string("python");
#endif
}

// Map helicity to which dvcsgen line to parse (matching your Python logic):
// pol=0  -> last line
// pol=+1 -> third from last
// pol=-1 -> second from last
// We will mirror exactly your historical usage: i = {0,2,1} counting from end.
static int dvcsgen_line_offset_from_end(Helicity h) {
    if (h == Helicity::Unpol) return 1; // last line
    if (h == Helicity::Plus)  return 3; // third from last
    return 2;                            // second from last
}

// Parse the last N lines and extract the numeric value at the end of the chosen line.
static double parse_dvcsgen_value(const std::string& all, int which_from_end) {
    // which_from_end: 1 means last line, 2 means second last, etc.
    std::vector<std::string> lines;
    {
        std::istringstream iss(all);
        std::string L;
        while (std::getline(iss, L)) {
            L = trim_copy(L);
            if (!L.empty()) lines.push_back(L);
        }
    }

    if ((int)lines.size() < which_from_end) return 0.0;

    const std::string& target = lines[lines.size() - which_from_end];

    std::istringstream ss(target);
    std::string tok;
    std::string last;
    while (ss >> tok) last = tok;

    try {
        return std::stod(last);
    } catch (...) {
        return 0.0;
    }
}

// Build dvcsgen command string.
static std::string build_dvcsgen_cmd(const std::string& dvcsgen_dir,
                                     double xB, double Q2, double t_pos, double phi_deg,
                                     double Ebeam, int bh_mode, bool globalfit) {
    // dvcsgen wants phi in radians.
    const double phi_rad = phi_deg * M_PI / 180.0;

    std::ostringstream cmd;
    cmd << dvcsgen_dir << "/dvcsgen"
        << " --beam " << Ebeam
        << " --x "    << xB << " " << xB
        << " --q2 "   << Q2 << " " << Q2
        << " --t "    << t_pos << " " << t_pos
        << " --bh "   << bh_mode
        << " --phi "  << phi_rad
        << " --ycol 0.0001";

    if (globalfit) cmd << " --globalfit";
    return cmd.str();
}

// Public: VGG total (BH+DVCS+INT) (your bh=3 case).
double vgg_xs(double xB, double Q2, double t_pos, double phi_deg,
              double Ebeam, Helicity helicity,
              const ModelPaths& paths, bool globalfit) {
    const std::string dvcsgen_dir = resolve_dvcsgen_dir(paths);

    std::ostringstream cmd;
    cmd << build_dvcsgen_cmd(dvcsgen_dir, xB, Q2, t_pos, phi_deg, Ebeam, /*bh_mode=*/3, globalfit)
        << " --gpd 101";

    std::ostringstream wrapped;
#ifndef _WIN32
    wrapped << "export PATH=" << shell_quote_posix(dvcsgen_dir) << ":$PATH; "
            << "export CLASDVCS_PDF=" << shell_quote_posix(dvcsgen_dir) << "; "
            << cmd.str();
#else
    wrapped << cmd.str();
#endif

    const std::string out = run_and_capture_stdout(wrapped.str(), /*require_ok=*/false);
    if (out.empty()) return 0.0;

    const int which = dvcsgen_line_offset_from_end(helicity);
    return parse_dvcsgen_value(out, which);
}

// Public: BH only.
double vgg_bh_only(double xB, double Q2, double t_pos, double phi_deg,
                   double Ebeam, const ModelPaths& paths, bool globalfit) {
    const std::string dvcsgen_dir = resolve_dvcsgen_dir(paths);

    std::ostringstream cmd;
    cmd << build_dvcsgen_cmd(dvcsgen_dir, xB, Q2, t_pos, phi_deg, Ebeam, /*bh_mode=*/1, globalfit);

    std::ostringstream wrapped;
#ifndef _WIN32
    wrapped << "export PATH=" << shell_quote_posix(dvcsgen_dir) << ":$PATH; "
            << "export CLASDVCS_PDF=" << shell_quote_posix(dvcsgen_dir) << "; "
            << cmd.str();
#else
    wrapped << cmd.str();
#endif

    const std::string out = run_and_capture_stdout(wrapped.str(), /*require_ok=*/false);
    if (out.empty()) return 0.0;

    // For BH-only you were taking the last line.
    return parse_dvcsgen_value(out, /*which_from_end=*/1);
}

// Public: KM15 via Python CLI.
// IMPORTANT: run with PYTHONPATH removed so ROOT module contamination cannot break numpy/gepard.
double km15_xs(double xB, double Q2, double t_pos, double phi_deg,
               double Ebeam, Helicity helicity,
               const ModelPaths& paths) {
    const std::string cli = resolve_km15_cli(paths);

    int hel_int = 0;
    if (helicity == Helicity::Plus)  hel_int = +1;
    if (helicity == Helicity::Minus) hel_int = -1;

#ifndef _WIN32
    const std::string py_km15 = resolve_km15_python();

    // Critical: drop PYTHONPATH for this subprocess only.  Also suppress stderr here:
    // km15_cli.py prints 0.0 on normal model errors, and repeated Python/env errors
    // otherwise flood long systematic jobs.
    std::ostringstream cmd;
    cmd << "env -u PYTHONPATH "
        << shell_quote_posix(py_km15) << " "
        << shell_quote_posix(cli)
        << " " << xB
        << " " << Q2
        << " " << t_pos
        << " " << phi_deg
        << " " << Ebeam
        << " " << hel_int
        << " XS";
#else
    std::ostringstream cmd;
    cmd << "python " << cli
        << " " << xB
        << " " << Q2
        << " " << t_pos
        << " " << phi_deg
        << " " << Ebeam
        << " " << hel_int
        << " XS";
#endif

    const std::string out = run_and_capture_stdout(cmd.str() + " 2>/dev/null", /*require_ok=*/false);
    if (out.empty()) return 0.0;

    try {
        return std::stod(trim_copy(out));
    } catch (...) {
        return 0.0;
    }
}