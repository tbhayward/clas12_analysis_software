#pragma once
#include <cmath>
#include <atomic>
#include <limits>
// Diagnostic-only Sebastian Table II bin assignment. Never changes exclusivity cuts.
inline std::atomic<bool> sebastian_use_tgamma_binning{false};
inline double sebastian_tgamma_abs(double ep, double eth, double eph,
                                   double gp, double gth, double gph,
                                   double beam_energy = 10.604) {
    constexpr double me = 0.00051099895;
    if (!(std::isfinite(ep) && std::isfinite(eth) && std::isfinite(eph) &&
          std::isfinite(gp) && std::isfinite(gth) && std::isfinite(gph) &&
          ep > 0 && gp > 0)) return std::numeric_limits<double>::quiet_NaN();
    const double ex=ep*std::sin(eth)*std::cos(eph), ey=ep*std::sin(eth)*std::sin(eph);
    const double ez=ep*std::cos(eth);
    const double gx=gp*std::sin(gth)*std::cos(gph), gy=gp*std::sin(gth)*std::sin(gph);
    const double gz=gp*std::cos(gth);
    const double q0=beam_energy-std::sqrt(ep*ep+me*me)-gp;
    const double qx=-ex-gx, qy=-ey-gy, qz=beam_energy-ez-gz;
    return std::fabs(q0*q0-qx*qx-qy*qy-qz*qz);
}
inline double sebastian_bin_t(double proton_t, double ep, double eth, double eph,
                               double gp, double gth, double gph) {
    return sebastian_use_tgamma_binning.load(std::memory_order_relaxed)
        ? sebastian_tgamma_abs(ep,eth,eph,gp,gth,gph) : std::fabs(proton_t);
}
