#!/usr/bin/env python3

import math

def calculate_dilution_error(nA, nC, nCH, nMT, nf, xA, xC, xCH, xHe, xf):
    # term1
    term1 = (5438.81 * nA * nf * xC**2 * xCH**2 * xf**2 * xHe**2
             * (-1.0 * nMT * xA + nA * xHe)**2
             * (1.0 * nMT * xC * xCH - 1.19655 * nCH * xC * xHe + 0.196547 * nC * xCH * xHe)**2)

    # term2
    term2 = (85044.9 * nA * nCH * xC**2 * xCH**2 * xf**2 * xHe**2
             * (-1.0 * nMT * xA + nA * xHe)**2
             * (1.0 * nMT * xC * xf - 0.302592 * nf * xC * xHe - 0.697408 * nC * xf * xHe)**2)

    # term3
    term3 = (47470.0 * nA * nC * xC**2 * xCH**2 * xf**2 * xHe**2
             * (-1.0 * nMT * xA + nA * xHe)**2
             * (1.0 * nMT * xCH * xf - 0.0665285 * nf * xCH * xHe - 0.933472 * nCH * xf * xHe)**2)

    # term4
    term4 = (1371.83 * nMT**2 * xA**2
             * (1.0 * nMT * xC * xCH * xf
                + (-1.97748 * nf * xC * xCH + 6.62306 * nCH * xC * xf - 5.64558 * nC * xCH * xf) * xHe)**2
             * (1.0 * nMT * xC * xCH * xf
                + (0.0136533 * nf * xC * xCH - 1.25054 * nCH * xC * xf + 0.236883 * nC * xCH * xf) * xHe)**2)

    # term5 numerator parts
    term5_num_p1 = 73.2426 * nMT**2 * xA * xC**2 * xCH**2 * xf**2
    term5_num_p2 = (nMT * xC * xCH * xf
                    * (2.0 * nf * xA * xC * xCH
                       + (-183.185 * nCH * xA * xC + 34.6998 * nC * xA * xCH + 1.42109e-14 * nA * xC * xCH) * xf))
    term5_num_p3 = ((-1.97748 * nf**2 * xA * xC**2 * xCH**2)
                    + nf * xC * xCH * (187.746 * nCH * xA * xC - 39.9547 * nC * xA * xCH - 145.836 * nA * xC * xCH) * xf
                    + (-606.623 * nCH**2 * xA * xC**2
                       + nCH * xC * (632.002 * nC * xA + 576.683 * nA * xC) * xCH
                       + nC * (-97.9502 * nC * xA - 430.847 * nA * xC) * xCH**2) * xf**2)

    term5_num = term5_num_p1 + term5_num_p2 * xHe + term5_num_p3 * xHe**2
    term5 = 0.255725 * nA * nMT * term5_num**2

    # denominator
    denom = (nA**3 * xHe**2
             * (73.2426 * nMT * xC * xCH * xf
                + 1.0 * nf * xC * xCH * xHe
                - 91.5925 * nCH * xC * xf * xHe
                + 17.3499 * nC * xCH * xf * xHe)**4)

    sigma_df = 27.3473 * math.sqrt((term1 + term2 + term3 + term4 + term5) / denom)
    return sigma_df

def calculate_dilution_and_error(nA, nC, nCH, nMT, nf, xA, xC, xCH, xHe, xf):
    dilution = ((27.3473 * (-1.0 * nMT * xA + nA * xHe)
                 * (-0.505693 * nMT * xC * xCH * xf
                    + (1.0 * nf * xC * xCH - 3.34924 * nCH * xC * xf + 2.85493 * nC * xCH * xf) * xHe))
                / (nA * xHe
                   * (73.2426 * nMT * xC * xCH * xf
                      + 1.0 * nf * xC * xCH * xHe
                      - 91.5925 * nCH * xC * xf * xHe
                      + 17.3499 * nC * xCH * xf * xHe)))
    error = calculate_dilution_error(nA, nC, nCH, nMT, nf, xA, xC, xCH, xHe, xf)
    return dilution, error

if __name__ == "__main__":
    # charge fractions (same for both cases)
    xA, xC, xCH, xHe, xf = 0.5737, 0.2057, 0.1855, 0.0291, 0.0061

    scenarios = {
        "8nA": {"nA": 143382290, "nC": 40762887, "nCH": 45881161, "nMT": 5914840, "nf": 126106},
        "4nA": {"nA": 149236000, "nC": 42425700, "nCH": 47754100, "nMT": 6156300, "nf": 131250},
    }

    for label, vals in scenarios.items():
        df, err = calculate_dilution_and_error(
            vals["nA"], vals["nC"], vals["nCH"], vals["nMT"], vals["nf"],
            xA, xC, xCH, xHe, xf
        )
        print(f"{label}: dilution = {df:.6f}, uncertainty = {err:.6f}")
    #endfor
#endif