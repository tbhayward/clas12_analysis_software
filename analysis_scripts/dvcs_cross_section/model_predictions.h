#ifndef MODEL_PREDICTIONS_H
#define MODEL_PREDICTIONS_H

#include <string>

// Simple enum for helicity selection where applicable.
// Use Helicity::Unpol for unpolarized, Helicity::Plus or Helicity::Minus for beam helicity dependent.
enum class Helicity {
    Unpol = 0,   // matches pol=0 
    Plus  = +1,  // matches pol=+1
    Minus = -1   // matches pol=-1
};

// Optional configuration object so callers can override paths without changing env vars.
struct ModelPaths {
    std::string dvcsgen_path;  // directory containing the dvcsgen executable
    std::string km15_cli;      // full path to km15_cli.py
    // If empty, the implementation will fall back to sensible defaults or env vars.
};

// ---------- Public API ----------
// All angles in degrees. t_pos is positive |t| (GeV^2). Beam energy in (GeV).

// VGG cross section (dvcsgen, "--bh 3" when you want BH+DVCS+interference "total").
// helicity controls which line of dvcsgen output we parse (as in your Python).
// Returns cross section value at the requested kinematics.
// Throws std::runtime_error on fatal errors; otherwise returns 0.0 if dvcsgen fails gracefully.
double vgg_xs(double xB, double Q2, double t_pos, double phi_deg,
              double Ebeam, Helicity helicity,
              const ModelPaths& paths = ModelPaths(), bool globalfit = false);

// Bethe-Heitler only (dvcsgen, "--bh 1").
double vgg_bh_only(double xB, double Q2, double t_pos, double phi_deg,
                   double Ebeam, const ModelPaths& paths = ModelPaths(),
                   bool globalfit = false);

// KM15 prediction via a thin Python CLI wrapper (gepard + th_KM15).
// Internally converts to phi_trento = pi - phi_rad to match your earlier usage.
double km15_xs(double xB, double Q2, double t_pos, double phi_deg,
               double Ebeam, Helicity helicity,
               const ModelPaths& paths = ModelPaths());

// ---------- Utilities ----------
// Set default search paths via environment variables (optional):
//   DVCSGEN_PATH : directory containing dvcsgen
//   KM15_CLI     : full path to km15_cli.py
//
// You can also pass ModelPaths when calling the functions. If a field is empty,
// the function will consult the environment, then fallback to a compiled-in default.
//
// Notes:
// - All functions expect t_pos >= 0 (they will send t = -|t_pos| to the model).
// - phi_deg is in degrees; conversions are handled internally.
// - For VGG helicity mapping: Unpol picks the "unpolarized total" line,
//   Plus/Minus pick the polarized-difference lines consistent with your Python selector.

#endif // MODEL_PREDICTIONS_H