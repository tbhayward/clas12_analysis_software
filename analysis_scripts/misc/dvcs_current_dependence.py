#!/usr/bin/env python3
#
# dvcs_current_dependence.py
#
# Hard-coded from the user-provided run-period summary:
#   - For each period and current: DVCS counts and accumulated charge (nC)
#   - Compute counts per nC with Poisson error bars
#   - Fit y = m*x + b (weighted least squares)
#   - Make two 2x3 canvases:
#       (1) counts/nC vs current with fits
#       (2) percent-of-intercept vs current, where intercept b is defined as 100%
#
# Updates requested by user:
#   - Chronological order: Sp18 Inb, Sp18 Out, Fa18 Inb, Fa18 Out, Sp19 Inb
#   - Drop Sp18 Inb points at 5 nA, 10 nA, and 75 nA
#   - Legends:
#       * remove "data (Poisson)" and "data / b (b -> 100)" text
#       * DO NOT use scientific notation; print all legend numbers to 5 decimals
#
# Output:
#   output/dvcs_counts_per_nc_vs_current.png
#   output/dvcs_percent_of_intercept_vs_current.png
#
# Run:
#   python3 dvcs_current_dependence.py
#

import os
import math
import numpy as np
import matplotlib.pyplot as plt


def sum_charge_expr(expr):
    """
    Sum a string expression like:
      "123.4+56.7+89.0"
    into a float.
    """
    if expr is None:
        return 0.0
    #endif

    s = str(expr).replace("\n", "").replace("\t", "").strip()
    if s == "":
        return 0.0
    #endif

    parts = s.split("+")
    total = 0.0
    for p in parts:
        pp = p.strip()
        if pp == "":
            continue
        #endif
        total += float(pp)
    #endfor

    return total
#endif


def weighted_linear_fit(x, y, sy):
    """
    Weighted least squares fit for y = m*x + b.

    Inputs:
      x, y, sy: numpy arrays
        - sy are 1-sigma uncertainties on y

    Returns dict:
      m, b, sm, sb, cov_mb, chi2, ndof
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sy = np.asarray(sy, dtype=float)

    if x.size < 2:
        return {
            "m": float("nan"),
            "b": float("nan"),
            "sm": float("nan"),
            "sb": float("nan"),
            "cov_mb": float("nan"),
            "chi2": float("nan"),
            "ndof": 0,
        }
    #endif

    w = 1.0 / (sy * sy)

    S = np.sum(w)
    Sx = np.sum(w * x)
    Sy = np.sum(w * y)
    Sxx = np.sum(w * x * x)
    Sxy = np.sum(w * x * y)

    D = S * Sxx - Sx * Sx
    if D == 0.0:
        return {
            "m": float("nan"),
            "b": float("nan"),
            "sm": float("nan"),
            "sb": float("nan"),
            "cov_mb": float("nan"),
            "chi2": float("nan"),
            "ndof": 0,
        }
    #endif

    m = (S * Sxy - Sx * Sy) / D
    b = (Sxx * Sy - Sx * Sxy) / D

    var_m = S / D
    var_b = Sxx / D
    cov_mb = -Sx / D

    sm = math.sqrt(var_m) if var_m >= 0.0 else float("nan")
    sb = math.sqrt(var_b) if var_b >= 0.0 else float("nan")

    yfit = m * x + b
    chi2 = np.sum(((y - yfit) / sy) ** 2)
    ndof = int(x.size - 2)

    return {
        "m": m,
        "b": b,
        "sm": sm,
        "sb": sb,
        "cov_mb": cov_mb,
        "chi2": chi2,
        "ndof": ndof,
    }
#endif


def period_points_from_dict(period_block, drop_currents=None):
    """
    period_block: dict mapping current (nA) -> {"counts": int, "charge": "a+b+c"}.

    drop_currents: optional set of currents (nA) to exclude.

    Returns:
      x (currents), y (counts/nC), sy (Poisson propagated), raw_counts, raw_charge
    """
    if drop_currents is None:
        drop_currents = set()
    #endif

    currents = sorted(period_block.keys())
    x = []
    y = []
    sy = []
    raw_counts = []
    raw_charge = []

    for I in currents:
        if I in drop_currents:
            continue
        #endif

        counts = float(period_block[I]["counts"])
        charge = sum_charge_expr(period_block[I]["charge"])

        if charge <= 0.0:
            raise RuntimeError(f"Charge <= 0 encountered for current {I} nA")
        #endif

        rate = counts / charge
        rate_err = math.sqrt(counts) / charge if counts > 0.0 else 0.0

        x.append(float(I))
        y.append(rate)
        sy.append(rate_err)
        raw_counts.append(counts)
        raw_charge.append(charge)
    #endfor

    return (
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(sy, dtype=float),
        np.asarray(raw_counts, dtype=float),
        np.asarray(raw_charge, dtype=float),
    )
#endif


def f5(val):
    """
    Format a float with 5 digits after the decimal, no scientific notation.
    """
    try:
        return f"{float(val):.5f}"
    except Exception:
        return "nan"
    #endif
#endif


def main():
    os.makedirs("output", exist_ok=True)

    data = {
        "Fa18 Inb": {
            5: {
                "counts": 2348,
                "charge": "18691.954+18474.5485+13923.970+13817.2078",
            },
            45: {
                "counts": 39075,
                "charge": "201346.38787+201166.16036+193579.079+194140.941+74896.2782+75326.12602+191423.3704+192464.25308+49888.8757+50064.795+66264.28332+66037.53247",
            },
            50: {
                "counts": 59234,
                "charge": "193667.2073+193616.37136+196779.37792+196386.098205+196626.8271+196294.92504+192402.17398+191979.34851+14640.9625+14810.92419+202062.1233+201417.8885",
            },
            55: {
                "counts": 48611,
                "charge": "52638.2290+52297.98364+3856.1777+3793.7708+188675.4283+188555.84970+185352.695+184093.51965+196134.55081+195204.362823+199060.8468+199012.79614+148373.1586+148479.6847+6656.71997+6607.21667",
            },
        },
        "Fa18 Out": {
            20: {"counts": 8811, "charge": "64156.3185+64410.056869"},
            40: {
                "counts": 144302,
                "charge": "143420.9199+143506.10577+145714.9450+146205.00158+35418.6917+35422.01965+75743.7664+76034.69329+133973.0374+133894.99798+144583.99597+144766.05413+144327.84197+144488.454528+143890.745208+144059.68832+137169.82696+137492.723388",
            },
            50: {
                "counts": 124850,
                "charge": "138224.19996+139038.25774+141704.8806+142387.018798+51007.15634+51297.464233+137540.88531+137841.00769+137045.91143+137018.58642+37793.01580+37947.76361+90860.10198+90944.84423+137721.85491+138138.37646+88172.906+88302.75451+8090.08282+8247.58453",
            },
            70: {"counts": 5197, "charge": "45040.08001+45273.50250"},
        },
        "Sp19 Inb": {
            5: {"counts": 867, "charge": "11384.66637+11332.34575"},
            10: {"counts": 3415, "charge": "67845.5772+67673.031265"},
            50: {
                "counts": 70165,
                "charge": "323344.4873+323692.674255+291200.93087+292563.48767+316595.28671+316518.92312+120891.25866+120788.63574+235832.3674+235687.911804+261436.456359+260771.4342041",
            },
        },
        "Sp18 Inb": {
            5: {"counts": 872, "charge": "2922.5921+3090.34080+2803.08867"},
            10: {"counts": 2504, "charge": "10381.05239+11344.23962"},
            35: {
                "counts": 22628,
                "charge": "143041.3698+180455.1287+99198.28764+215819.24902+253604.037+256674.04077",
            },
            45: {"counts": 3649, "charge": "247575.2544"},
            50: {
                "counts": 29985,
                "charge": "86823.20456+53875.207720000006+462365.5763500001+155942.65868+30313.130000000005+190978.2245096+57414.873250000004+299534.5488+15610.063702499998+336045.59401999996",
            },
            55: {
                "counts": 25812,
                "charge": "494237.94639999996+492764.5799800001+68830.43074000001+253110.41491+236823.99630479995+61038.981349999995",
            },
            75: {
                "counts": 2406,
                "charge": "7457.454999999996+6046.243000000002+10627.161851500001",
            },
        },
        "Sp18 Out": {
            25: {"counts": 6089, "charge": "107432.82072000002"},
            30: {
                "counts": 50486,
                "charge": "113278.01032+11124.354809999999+93134.14310000002+128407.92089999998+34088.218689999994+55631.41555+87828.17699999998+125198.45608000002+114118.60385999999+95226.42483999998+107480.42574000002",
            },
            45: {
                "counts": 92032,
                "charge": "283283.2660000002+53094.88620000001+296750.10305000003+292795.32966000005+146204.13000000003+294288.15236000007+259937.71113000007+129771.95737+290971.9334",
            },
        },
    }

    period_order = ["Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"]

    drop_currents_by_period = {
        "Sp18 Inb": {5, 10, 75},
    }

    fig1, axs1 = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
    axs1 = axs1.flatten()

    fit_results = {}
    xfit = np.linspace(0.0, 100.0, 300)

    print("")
    print("=== Canvas 1 fits: y = m*x + b for counts_per_nC vs current ===")

    for i, period in enumerate(period_order):
        ax = axs1[i]

        drop = drop_currents_by_period.get(period, set())
        x, y, sy, raw_counts, raw_charge = period_points_from_dict(data[period], drop_currents=drop)

        fr = weighted_linear_fit(x, y, sy)
        fit_results[period] = fr

        m = fr["m"]
        b = fr["b"]
        sm = fr["sm"]
        sb = fr["sb"]
        chi2 = fr["chi2"]
        ndof = fr["ndof"]

        print(f"{period}: m = {m:.10f} +/- {sm:.10f}, b = {b:.10f} +/- {sb:.10f}, chi2/ndof = {chi2:.2f}/{ndof}")

        eb = ax.errorbar(
            x, y, yerr=sy,
            fmt="o",
            capsize=3,
            label=f"m={f5(m)} +/- {f5(sm)}\n b={f5(b)} +/- {f5(sb)}"
        )
        color = eb[0].get_color()
        ax.plot(xfit, m * xfit + b, color=color)

        ax.set_title(period)
        ax.set_xlim(0.0, 100.0)
        ax.set_xlabel("Beam current (nA)")
        ax.set_ylabel("DVCS counts per charge (1/nC)")
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=True)
    #endfor

    axc = axs1[5]
    axc.set_title("All periods (overlay)")
    axc.set_xlim(0.0, 100.0)
    axc.set_xlabel("Beam current (nA)")
    axc.set_ylabel("DVCS counts per charge (1/nC)")
    axc.grid(True, alpha=0.3)

    for period in period_order:
        drop = drop_currents_by_period.get(period, set())
        x, y, sy, raw_counts, raw_charge = period_points_from_dict(data[period], drop_currents=drop)

        fr = fit_results[period]
        m = fr["m"]
        b = fr["b"]
        sm = fr["sm"]
        sb = fr["sb"]

        eb = axc.errorbar(
            x, y, yerr=sy,
            fmt="o",
            capsize=3,
            label=f"{period}: m={f5(m)} +/- {f5(sm)}, b={f5(b)} +/- {f5(sb)}"
        )
        color = eb[0].get_color()
        axc.plot(xfit, m * xfit + b, color=color)
    #endfor

    axc.legend(frameon=True, fontsize=9)

    out1 = "output/dvcs_counts_per_nc_vs_current.png"
    fig1.savefig(out1, dpi=200)
    print("")
    print(f"[saved] {out1}")

    fig2, axs2 = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
    axs2 = axs2.flatten()

    print("")
    print("=== Canvas 2 derived slopes: percent(x) = 100*(m*x+b)/b ===")

    for i, period in enumerate(period_order):
        ax = axs2[i]

        drop = drop_currents_by_period.get(period, set())
        x, y, sy, raw_counts, raw_charge = period_points_from_dict(data[period], drop_currents=drop)

        fr = fit_results[period]
        m = fr["m"]
        b = fr["b"]
        sm = fr["sm"]
        sb = fr["sb"]
        cov_mb = fr["cov_mb"]

        pct = 100.0 * (y / b)
        pct_err = 100.0 * np.sqrt((sy / b) ** 2 + ((y * sb) / (b * b)) ** 2)
        pct_fit = 100.0 * ((m * xfit + b) / b)

        slope_pct = 100.0 * (m / b)

        dp_dm = 100.0 / b
        dp_db = -100.0 * m / (b * b)
        var_p = dp_dm * dp_dm * (sm * sm) + dp_db * dp_db * (sb * sb) + 2.0 * dp_dm * dp_db * cov_mb
        slope_pct_err = math.sqrt(var_p) if var_p >= 0.0 else float("nan")

        print(f"{period}: b = {b:.10f} +/- {sb:.10f}, slope = {slope_pct:.10f} +/- {slope_pct_err:.10f} (% per nA)")

        eb = ax.errorbar(
            x, pct, yerr=pct_err,
            fmt="o",
            capsize=3,
            label=f"b={f5(b)} +/- {f5(sb)}\n slope={f5(slope_pct)} +/- {f5(slope_pct_err)} (%/nA)"
        )
        color = eb[0].get_color()
        ax.plot(xfit, pct_fit, color=color)

        ax.set_title(period)
        ax.set_xlim(0.0, 100.0)
        ax.set_xlabel("Beam current (nA)")
        ax.set_ylabel("Percent of intercept b (%)")
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=True)
    #endfor

    axc2 = axs2[5]
    axc2.set_title("All periods (overlay)")
    axc2.set_xlim(0.0, 100.0)
    axc2.set_xlabel("Beam current (nA)")
    axc2.set_ylabel("Percent of intercept b (%)")
    axc2.grid(True, alpha=0.3)

    for period in period_order:
        drop = drop_currents_by_period.get(period, set())
        x, y, sy, raw_counts, raw_charge = period_points_from_dict(data[period], drop_currents=drop)

        fr = fit_results[period]
        m = fr["m"]
        b = fr["b"]
        sm = fr["sm"]
        sb = fr["sb"]
        cov_mb = fr["cov_mb"]

        pct = 100.0 * (y / b)
        pct_err = 100.0 * np.sqrt((sy / b) ** 2 + ((y * sb) / (b * b)) ** 2)
        pct_fit = 100.0 * ((m * xfit + b) / b)

        slope_pct = 100.0 * (m / b)
        dp_dm = 100.0 / b
        dp_db = -100.0 * m / (b * b)
        var_p = dp_dm * dp_dm * (sm * sm) + dp_db * dp_db * (sb * sb) + 2.0 * dp_dm * dp_db * cov_mb
        slope_pct_err = math.sqrt(var_p) if var_p >= 0.0 else float("nan")

        eb = axc2.errorbar(
            x, pct, yerr=pct_err,
            fmt="o",
            capsize=3,
            label=f"{period}: slope={f5(slope_pct)} +/- {f5(slope_pct_err)} (%/nA)"
        )
        color = eb[0].get_color()
        axc2.plot(xfit, pct_fit, color=color)
    #endfor

    axc2.legend(frameon=True, fontsize=9)

    out2 = "output/dvcs_percent_of_intercept_vs_current.png"
    fig2.savefig(out2, dpi=200)
    print("")
    print(f"[saved] {out2}")
    print("")
#endif


if __name__ == "__main__":
    main()
#endif