#!/usr/bin/env python3
"""
calculate_dilution_factor.py

Module to compute dilution factor Df(xB) and its uncertainty
for each xB bin for RGC Su22 and RGC Fa22 datasets.
Saves CSV and PDF outputs, with verbose logging for progress tracking.
Filters events by y < 0.75 before counting.
"""

import numpy as np
import pandas as pd
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib.pyplot as plt
import os

# Updated charge fractions for each run period based on measured totals:
# RGC Su22: NH3=72.39%, C=7.14%, CH2=3.73%, He=7.51%, ET=9.23%
# RGC Fa22: NH3=56.05%, C=20.42%, CH2=18.66%, He=4.29%, ET=0.58%
CHARGE_FRAC = {
    "RGC_Su22": {"xNH3": 0.7239, "xC": 0.0714, "xCH2": 0.0373, "xHe": 0.0523, "xET": 0.1151},
    "RGC_Fa22": {"xNH3": 0.5690, "xC": 0.2072, "xCH2": 0.1894, "xHe": 0.0285, "xET": 0.0059},
    "RGC_Sp23": {"xNH3": 0.5629, "xC": 0.1331, "xCH2": 0.1517, "xHe": 0.0926, "xET": 0.0597},
}

# y-cut threshold
Y_CUT = 0.75


def calculate_dilution_error(nNH3, nC, nCH2, nHe, nET, xNH3, xC, xCH2, xHe, xET):
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
            (1.0 * nHe * xC * xCH2 * xET + (-1.97748 * nET * xC * xCH2 + 6.62306 * nCH2 * xC * xET - 5.64558 * nC * xCH2 * xET) * xHe)**2 * \
            (1.0 * nHe * xC * xCH2 * xET + (0.0136533 * nET * xC * xCH2 - 1.25054 * nCH2 * xC * xET + 0.236883 * nC * xCH2 * xET) * xHe)**2
    t5_p1 = 73.2426 * nHe**2 * xNH3 * xC**2 * xCH2**2 * xET**2
    t5_p2 = nHe * xC * xCH2 * xET * \
            (2.0 * nET * xNH3 * xC * xCH2 + (-183.185 * nCH2 * xNH3 * xC + 34.6998 * nC * xNH3 * xCH2 + 1.42109e-14 * nNH3 * xC * xCH2) * xET)
    t5_p3 = (-1.97748 * nET**2 * xNH3 * xC**2 * xCH2**2) + \
             nET * xC * xCH2 * (187.746 * nCH2 * xNH3 * xC - 39.9547 * nC * xNH3 * xCH2 - 145.836 * nNH3 * xC * xCH2) * xET + \
             (-606.623 * nCH2**2 * xNH3 * xC**2 + nCH2 * xC * (632.002 * nC * xNH3 + 576.683 * nNH3 * xC) * xCH2 + \
             nC * (-97.9502 * nC * xNH3 - 430.847 * nNH3 * xC) * xCH2**2) * xET**2
    term5 = 0.255725 * nNH3 * nHe * (t5_p1 + t5_p2 * xHe + t5_p3 * xHe**2)**2
    denom = nNH3**3 * xHe**2 * \
            (73.2426 * nHe * xC * xCH2 * xET + 1.0 * nET * xC * xCH2 * xHe - 91.5925 * nCH2 * xC * xET * xHe + 17.3499 * nC * xCH2 * xET * xHe)**4
    total = term1 + term2 + term3 + term4 + term5
    print(f"[Debug Error] sum_terms={total:.4e}, denom={denom:.4e}")
    var = total/denom if denom else 0.0
    if var < 0:
        print(f"[Warning] Negative var={var:.4e}, clamping to 0.")
        var = 0.0
    sigma_df = 27.3473 * math.sqrt(var)
    print(f"[Error] sigma_Df={sigma_df:.4e}")
    return sigma_df


def calculate_dilution_and_error(nNH3, nC, nCH2, nHe, nET, xNH3, xC, xCH2, xHe, xET):
    print(f"[Dilution] Df counts: NH3={nNH3}, C={nC}, CH2={nCH2}, He={nHe}, ET={nET}")
    num = (27.3473 * (-1.0 * nHe * xNH3 + nNH3 * xHe) *
           (-0.505693 * nHe * xC * xCH2 * xET
              + (1.0 * nET * xC * xCH2 - 3.34924 * nCH2 * xC * xET + 2.85493 * nC * xCH2 * xET) * xHe))
    denom = (nNH3 * xHe
             * (73.2426 * nHe * xC * xCH2 * xET
                + 1.0 * nET * xC * xCH2 * xHe
                - 91.5925 * nCH2 * xC * xET * xHe
                + 17.3499 * nC * xCH2 * xET * xHe))
    print(f"[Debug Dilution] num={num:.4e}, denom={denom:.4e}")
    df = num/denom if denom else 0.0
    print(f"[Dilution] Df={df:.4e}")
    err = calculate_dilution_error(nNH3, nC, nCH2, nHe, nET, xNH3, xC, xCH2, xHe, xET)
    print(f"[Dilution] Df_err={err:.4e}\n")
    return df, err


def _compute_period(period, trees, xB_bins):
    print(f"[Period {period}] Starting with y<{Y_CUT}")
    data_x, data_y = {}, {}
    for t, tr in trees.items():
        data_x[t] = tr["x"].array(library="np")
        data_y[t] = tr["y"].array(library="np")
    nbins = len(xB_bins)-1
    rows = []
    fr = CHARGE_FRAC[period]
    for i in range(nbins):
        idx = {t: np.digitize(data_x[t][data_y[t]<Y_CUT], xB_bins)-1 for t in data_x}
        x_mean = np.mean(data_x["NH3"][data_y["NH3"]<Y_CUT][idx["NH3"]==i]) if np.any(idx["NH3"]==i) else 0.0
        counts = {t: np.sum(idx[t]==i) for t in idx}
        print(f"[Period {period}] Bin{i+1}/{nbins} x_mean={x_mean:.4f} counts={counts}")
        Df, Df_s = calculate_dilution_and_error(
            counts["NH3"],counts["C"],counts["CH2"],counts["He"],counts["ET"],
            fr["xNH3"],fr["xC"],fr["xCH2"],fr["xHe"],fr["xET"]
        )
        rows.append({"x_mean": x_mean, "Df": Df, "Df_sigma": Df_s})
    df = pd.DataFrame(rows)
    df["period"] = period
    print(f"[Period {period}] Done.\n")
    return df


def calculate_and_save(trees, xB_bins):
    print("[Main] Begin dilution calculation")
    os.makedirs("output", exist_ok=True)
    per = ["RGC_Su22", "RGC_Fa22"]
    with ProcessPoolExecutor() as ex:
        futs = {ex.submit(_compute_period, p, trees[p], xB_bins): p for p in per}
        res = {futs[f]: f.result() for f in as_completed(futs)}
    df_su, df_fa = res["RGC_Su22"], res["RGC_Fa22"]
    wsu = 1/df_su["Df_sigma"]**2; wfa = 1/df_fa["Df_sigma"]**2
    df_avg = pd.DataFrame({
        "x_mean": df_su["x_mean"],
        "Df": (df_su["Df"]*wsu + df_fa["Df"]*wfa)/(wsu+wfa),
        "Df_sigma": np.sqrt(1/(wsu+wfa)),
        "period": "weighted"
    })
    df_all = pd.concat([df_su, df_fa, df_avg], ignore_index=True)
    csv = "output/dilution_factor.csv"; df_all.to_csv(csv, index=False); print(f"Saved CSV {csv}")
    plt.errorbar(df_su["x_mean"], df_su["Df"], yerr=df_su["Df_sigma"], fmt="o", label="Su22")
    plt.errorbar(df_fa["x_mean"], df_fa["Df"], yerr=df_fa["Df_sigma"], fmt="o", label="Fa22")
    plt.xlabel("$x_B$"); plt.ylabel("$D_f$"); plt.legend(loc="upper right"); plt.xlim(0,0.8); plt.ylim(0.1,0.3); plt.tight_layout()
    pdf = "output/dilution_factor.pdf"; plt.savefig(pdf); plt.close(); print(f"Saved PDF {pdf}\nDone.")

# ------------------------------------------------------------------
# Temporary routine using manual dilution values provided by user
# ------------------------------------------------------------------

def calculate_dilution_factor_temp(trees, xB_bins):
    """
    Use manually provided Df values for Su22, Fa22, Sp23,
    compute combined weighted average, uncertainty, population std,
    and range, then save CSV and plot.
    Mean xB per bin is calculated from the NH3 tree with y<Y_CUT.
    """
    # Load NH3 x and y for all three periods
    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    data_x = {p: trees[p]["NH3"]["x"].array(library="np") for p in periods}
    data_y = {p: trees[p]["NH3"]["y"].array(library="np") for p in periods}
    nbins = len(xB_bins) - 1

    # NEW Manual values for each period and bin
    su22_df  = np.array([0.166709, 0.176158, 0.191176, 0.206910, 0.221360, 0.230502, 0.220180])
    su22_err = np.array([0.000938862, 0.000491330, 0.000561747, 0.000733391, 0.001159140, 0.002228200, 0.005815470])

    fa22_df  = np.array([0.139084, 0.139782, 0.155144, 0.171282, 0.184396, 0.194965, 0.191098])
    fa22_err = np.array([0.000345853, 0.000188008, 0.000215951, 0.000283221, 0.000452165, 0.000864904, 0.002219650])

    sp23_df  = np.array([0.154322, 0.162926, 0.175352, 0.187502, 0.199844, 0.214446, 0.218649])
    sp23_err = np.array([0.000735081, 0.000389901, 0.000453347, 0.000602108, 0.000960937, 0.001828120, 0.004582990])

    # Compute mean xB per bin using Su22 NH3 and the y<Y_CUT cut
    mask = data_x["RGC_Su22"][data_y["RGC_Su22"] < Y_CUT]
    bidx = np.digitize(mask, xB_bins) - 1
    x_mean = np.array([
        np.mean(mask[bidx == i]) if np.any(bidx == i) else 0.0
        for i in range(nbins)
    ])

    # Weighted average across three periods
    w_su = 1.0 / su22_err**2
    w_fa = 1.0 / fa22_err**2
    w_sp = 1.0 / sp23_err**2
    w_sum = w_su + w_fa + w_sp
    df_avg = (su22_df*w_su + fa22_df*w_fa + sp23_df*w_sp) / w_sum
    err_avg = np.sqrt(1.0 / w_sum)

    # Population standard deviation across the three
    df_std = np.sqrt(
        ((su22_df - df_avg)**2 +
         (fa22_df - df_avg)**2 +
         (sp23_df - df_avg)**2) / 3.0
    )

    # Range (max - min) across three periods for each bin
    df_range = np.max([su22_df, fa22_df, sp23_df], axis=0) - \
               np.min([su22_df, fa22_df, sp23_df], axis=0)

    # Assemble DataFrame
    df_temp = pd.DataFrame({
        'x_mean':   x_mean,
        'Df_Su22':  su22_df,
        'Err_Su22': su22_err,
        'Df_Fa22':  fa22_df,
        'Err_Fa22': fa22_err,
        'Df_Sp23':  sp23_df,
        'Err_Sp23': sp23_err,
        'Df_avg':   df_avg,
        'Err_avg':  err_avg,
        'Df_std':   df_std,
        'Df_range': df_range,
    })

    # Save CSV
    os.makedirs('output', exist_ok=True)
    csv_path = 'output/dilution_factor.csv'
    df_temp.to_csv(csv_path, index=False)
    print(f"[Temp] Saved manual dilution factors CSV to {csv_path}")

    # Plot all three
    plt.errorbar(x_mean, su22_df,  yerr=su22_err, fmt='o', label='Su22')
    plt.errorbar(x_mean, fa22_df,  yerr=fa22_err, fmt='s', label='Fa22')
    plt.errorbar(x_mean, sp23_df, yerr=sp23_err, fmt='^', label='Sp23')
    plt.xlabel('$x_{B}$')
    plt.ylabel('$D_{f}$')
    plt.xlim(0, 0.8)
    plt.ylim(0.1, 0.3)
    plt.legend(loc='upper right')
    plt.tight_layout()
    pdf_path = 'output/dilution_factor.pdf'
    plt.savefig(pdf_path)
    plt.close()
    print(f"[Temp] Saved manual dilution factors plot to {pdf_path}")
