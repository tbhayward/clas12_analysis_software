#pragma once
// Lee pass-1 note (final April 2025), Sec. 5.4 and embedded paper Eq. (10).
// Frad = <radiative>/<Born>; Fbin = <Born>/Born(reference).
namespace dvcs_corrections {
inline double cross_section(double yield, double luminosity_volume,
                            double frad, double fbin) {
    return yield / (luminosity_volume * frad * fbin);
}
}
