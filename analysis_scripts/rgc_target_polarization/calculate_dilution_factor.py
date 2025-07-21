#!/usr/bin/env python3
"""
calculate_dilution_factor.py

Module to compute dilution factor Df(xB) and its uncertainty
for each xB bin for RGC Su22, Fa22 and Sp23 datasets.
Saves CSV and PDF outputs, with verbose logging for progress tracking.
Filters events by y < 0.75 before counting.
Supports an optional use_cached flag to pull x_mean from the CSV.
"""

import numpy as np
import pandas as pd
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib.pyplot as plt
import os

# Updated charge fractions for each run period based on measured totals:
CHARGE_FRAC = {
    "RGC_Su22": {"xNH3": 0.7452, "xC": 0.0735, "xCH2": 0.0383, "xHe": 0.0245, "xET": 0.1184},
    "RGC_Fa22": {"xNH3": 0.5737, "xC": 0.2057, "xCH2": 0.1855, "xHe": 0.0291, "xET": 0.0061},
    "RGC_Sp23": {"xNH3": 0.5632, "xC": 0.1331, "xCH2": 0.1514, "xHe": 0.0926, "xET": 0.0597},
}

# y-cut threshold
Y_CUT = 0.75


def calculate_dilution_error(nNH3, nC, nCH2, nHe, nET,
                             xNH3, xC, xCH2, xHe, xET):
    """
    Compute sigma_Df with debug output of each term.
    """
    a = -1.0 * nHe * xNH3 + nNH3 * xHe

    term1 = 5438.81 * nNH3 * nET * xC**2 * xCH2**2 * xET**2 * xHe**2 * a**2 * \
            (1.0 * nHe * xC * xCH2 - 1.19655 * nCH2 * xC * xHe + 0.196547 * nC * xCH2 * xHe)**2

    term2 = 85044.9 * nNH3 * nCH2 * xC**2 * xCH2**2 * xET**2 * xHe**2 * a**2 * \
            (1.0 * nHe * xC * xET - 0.302592 * nET * xC * xHe - 0.697408 * nC * xET * xHe)**2

    term3 = 47470.0 * nNH3 * nC * xC**2 * xCH2**2 * xET**2 * xHe**2 * a**2 * \
            (1.0 * nHe * xCH2 * xET - 0.0665285 * nET * xCH2 * xHe - 0.933472 * nCH2 * xET * xHe)**2

    term4 = 1371.83 * nHe**2 * xNH3**2 * \
            (1.0 * nHe * xC * xCH2 * xET + (-1.97748 * nET * xC * xCH2 +
             6.62306 * nCH2 * xC * xET - 5.64558 * nC * xCH2 * xET) * xHe)**2 * \
            (1.0 * nHe * xC * xCH2 * xET + (0.0136533 * nET * xC * xCH2 -
             1.25054 * nCH2 * xC * xET + 0.236883 * nC * xCH2 * xET) * xHe)**2

    t5_p1 = 73.2426 * nHe**2 * xNH3 * xC**2 * xCH2**2 * xET**2

    t5_p2 = nHe * xC * xCH2 * xET * \
            (2.0 * nET * xNH3 * xC * xCH2 +
             (-183.185 * nCH2 * xNH3 * xC +
              34.6998 * nC * xNH3 * xCH2 +
              1.42109e-14 * nNH3 * xC * xCH2) * xET)

    t5_p3 = (-1.97748 * nET**2 * xNH3 * xC**2 * xCH2**2) + \
             nET * xC * xCH2 * (187.746 * nCH2 * xNH3 * xC -
                                39.9547 * nC * xNH3 * xCH2 -
                                145.836 * nNH3 * xC * xCH2) * xET + \
             (-606.623 * nCH2**2 * xNH3 * xC**2 +
              nCH2 * xC * (632.002 * nC * xNH3 +
                           576.683 * nNH3 * xC) * xCH2 +
              nC * (-97.9502 * nC * xNH3 -
                     430.847 * nNH3 * xC) * xCH2**2) * xET**2

    term5 = 0.255725 * nNH3 * nHe * (t5_p1 + t5_p2 * xHe + t5_p3 * xHe**2)**2

    denom = nNH3**3 * xHe**2 * \
            (73.2426 * nHe * xC * xCH2 * xET +
             1.0 * nET * xC * xCH2 * xHe -
             91.5925 * nCH2 * xC * xET * xHe +
             17.3499 * nC * xCH2 * xET * xHe)**4

    total = term1 + term2 + term3 + term4 + term5
    var = total / denom if denom else 0.0
    if var < 0:
        var = 0.0
    sigma_df = 27.3473 * math.sqrt(var)
    return sigma_df


def calculate_dilution_and_error(nNH3, nC, nCH2, nHe, nET,
                                 xNH3, xC, xCH2, xHe, xET):
    """
    Compute the dilution factor Df and its error for given counts and charge fractions.
    """
    num = (
        27.3473
        * (-1.0 * nHe * xNH3 + nNH3 * xHe)
        * (
            -0.505693 * nHe * xC * xCH2 * xET
            + (1.0 * nET * xC * xCH2
               - 3.34924 * nCH2 * xC * xET
               + 2.85493 * nC * xCH2 * xET) * xHe
        )
    )

    denom = (
        nNH3
        * xHe
        * (
            73.2426 * nHe * xC * xCH2 * xET
            + 1.0 * nET * xC * xCH2 * xHe
            - 91.5925 * nCH2 * xC * xET * xHe
            + 17.3499 * nC * xCH2 * xET * xHe
        )
    )

    df = num / denom if denom else 0.0
    err = calculate_dilution_error(nNH3, nC, nCH2, nHe, nET,
                                   xNH3, xC, xCH2, xHe, xET)
    return df, err


def _compute_period(period, trees, xB_bins):
    """
    Internal: compute Df and error for one period by scanning the ROOT trees.
    """
    print(f"[Period {period}] Starting with y<{Y_CUT}")
    data_x, data_y = {}, {}
    for tgt, tr in trees.items():
        data_x[tgt] = tr["x"].array(library="np")
        data_y[tgt] = tr["y"].array(library="np")

    nbins = len(xB_bins) - 1
    rows = []
    fr = CHARGE_FRAC[period]

    for i in range(nbins):
        idx = {tgt: np.digitize(data_x[tgt][data_y[tgt] < Y_CUT], xB_bins) - 1
               for tgt in data_x}

        x_mean = (
            np.mean(data_x["NH3"][data_y["NH3"] < Y_CUT][idx["NH3"] == i])
            if np.any(idx["NH3"] == i)
            else 0.0
        )
        counts = {tgt: np.sum(idx[tgt] == i) for tgt in idx}

        Df, Df_s = calculate_dilution_and_error(
            counts["NH3"], counts["C"], counts["CH2"], counts["He"], counts["ET"],
            fr["xNH3"], fr["xC"], fr["xCH2"], fr["xHe"], fr["xET"]
        )
        rows.append({"x_mean": x_mean, "Df": Df, "Df_sigma": Df_s})

    df = pd.DataFrame(rows)
    df["period"] = period
    print(f"[Period {period}] Done.\n")
    return df


def calculate_and_save(trees, xB_bins):
    """
    Full calculation: uses ProcessPoolExecutor to run _compute_period
    for Su22 and Fa22, then forms a weighted average.
    """
    print("[Main] Begin dilution calculation")
    os.makedirs("output", exist_ok=True)

    periods = ["RGC_Su22", "RGC_Fa22"]
    with ProcessPoolExecutor() as ex:
        futures = {ex.submit(_compute_period, p, trees[p], xB_bins): p
                   for p in periods}
        results = {period: fut.result() for fut, period in futures.items()}

    df_su = results["RGC_Su22"]
    df_fa = results["RGC_Fa22"]

    wsu = 1.0 / df_su["Df_sigma"]**2
    wfa = 1.0 / df_fa["Df_sigma"]**2

    df_avg = pd.DataFrame({
        "x_mean": df_su["x_mean"],
        "Df":      (df_su["Df"] * wsu + df_fa["Df"] * wfa) / (wsu + wfa),
        "Df_sigma": np.sqrt(1.0 / (wsu + wfa)),
        "period":  "weighted",
    })

    df_all = pd.concat([df_su, df_fa, df_avg], ignore_index=True)
    df_all.to_csv("output/dilution_factor.csv", index=False)
    print("[Main] Saved CSV output/dilution_factor.csv")

    plt.errorbar(df_su["x_mean"], df_su["Df"],  yerr=df_su["Df_sigma"], fmt="o", label="Su22")
    plt.errorbar(df_fa["x_mean"], df_fa["Df"],  yerr=df_fa["Df_sigma"], fmt="o", label="Fa22")
    plt.xlabel("$x_B$"); plt.ylabel("$D_f$")
    plt.xlim(0, 0.8); plt.ylim(0.1, 0.3)
    plt.legend(loc="upper right"); plt.tight_layout()
    plt.savefig("output/dilution_factor.pdf"); plt.close()
    print("[Main] Saved PDF output/dilution_factor.pdf\nDone.")


def calculate_dilution_factor_temp(trees, xB_bins, use_cached=False):
    """
    Use manually provided Df values for Su22, Fa22, Sp23,
    compute combined weighted average, uncertainty, pop std and range,
    then save CSV and plot.
    If use_cached==True, pull x_mean from existing CSV instead of recalculating.
    """

    su22_df  = np.array([0.171112, 0.180644, 0.195561, 0.211246, 0.225397, 0.234519, 0.224291])
    su22_err = np.array([0.000927088, 0.000484202, 0.000552962, 0.000720490, 0.001138660, 0.002188790, 0.005721290])

    fa22_df  = np.array([0.139084, 0.139782, 0.155144, 0.171282, 0.184396, 0.194965, 0.191098])
    fa22_err = np.array([0.000345853, 0.000188008, 0.000215951, 0.000283221, 0.000452165, 0.000864904, 0.002219650])

    sp23_df  = np.array([0.156379, 0.164955, 0.177343, 0.189456, 0.201761, 0.216320, 0.220511])
    sp23_err = np.array([0.000732864, 0.000388726, 0.000451983, 0.000600302, 0.000958067, 0.001822700, 0.004569400])

    nbins = len(xB_bins) - 1

    # decide whether to recalc x_mean or load from cache
    if use_cached and os.path.exists("output/dilution_factor.csv"):
        df_cached = pd.read_csv("output/dilution_factor.csv")
        x_mean    = df_cached["x_mean"].to_numpy()
        print("[Temp] Using cached x_mean from CSV")
    else:
        # re‐calculate just once from Su22 NH3
        arr_x = trees["RGC_Su22"]["NH3"]["x"].array(library="np")
        arr_y = trees["RGC_Su22"]["NH3"]["y"].array(library="np")
        mask  = arr_x[arr_y < Y_CUT]
        idx   = np.digitize(mask, xB_bins) - 1
        x_mean = np.array([
            mask[idx == i].mean() if np.any(idx == i) else 0.0
            for i in range(nbins)
        ])

    # weighted average
    w_su = 1.0 / su22_err**2
    w_fa = 1.0 / fa22_err**2
    w_sp = 1.0 / sp23_err**2
    w_tot = w_su + w_fa + w_sp

    df_avg  = (su22_df*w_su + fa22_df*w_fa + sp23_df*w_sp) / w_tot
    err_avg = np.sqrt(1.0 / w_tot)

    # population standard deviation
    df_std = np.sqrt(
        ((su22_df - df_avg)**2 +
         (fa22_df - df_avg)**2 +
         (sp23_df - df_avg)**2) / 3.0
    )

    # range
    df_range = np.max([su22_df, fa22_df, sp23_df], axis=0) - \
               np.min([su22_df, fa22_df, sp23_df], axis=0)

    df_temp = pd.DataFrame({
        'x_mean':   x_mean,
        'Df_Su22':  su22_df,  'Err_Su22': su22_err,
        'Df_Fa22':  fa22_df,  'Err_Fa22': fa22_err,
        'Df_Sp23':  sp23_df,  'Err_Sp23': sp23_err,
        'Df_avg':   df_avg,   'Err_avg':  err_avg,
        'Df_std':   df_std,   'Df_range': df_range,
    })

    os.makedirs("output", exist_ok=True)
    df_temp.to_csv("output/dilution_factor.csv", index=False)
    print("[Temp] Saved manual dilution factors CSV to output/dilution_factor.csv")

    # plot all three
    plt.errorbar(x_mean, su22_df,  yerr=su22_err, fmt='o', label='Su22')
    plt.errorbar(x_mean, fa22_df,  yerr=fa22_err, fmt='s', label='Fa22')
    plt.errorbar(x_mean, sp23_df, yerr=sp23_err, fmt='^', label='Sp23')
    plt.xlabel('$x_{B}$'); plt.ylabel('$D_{f}$')
    plt.xlim(0, 0.8); plt.ylim(0.1, 0.3)
    plt.legend(loc='upper right'); plt.tight_layout()
    plt.savefig("output/dilution_factor.pdf"); plt.close()
    print("[Temp] Saved manual dilution factors plot to output/dilution_factor.pdf")