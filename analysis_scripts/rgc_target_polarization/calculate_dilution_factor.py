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
# RGC Fa22: NH3=58.49%, C=21.30%, CH2=19.47%, He=0.13%, ET=0.61%
CHARGE_FRAC = {
    "RGC_Su22": {"xNH3": 0.7239, "xC": 0.0714, "xCH2": 0.0373, "xHe": 0.0751, "xET": 0.0923},
    "RGC_Fa22": {"xNH3": 0.5849, "xC": 0.2130, "xCH2": 0.1947, "xHe": 0.0013, "xET": 0.0061},
}

# y-cut threshold
Y_CUT = 0.75

def calculate_dilution_error(nNH3, nC, nCH2, nHe, nET, xNH3, xC, xCH2, xHe, xET):
    """
    Compute sigma_Df with debug output of each term.
    """
    # a = -1.0 * nHe * xNH3 + nNH3 * xHe
    # term1 = 5438.81 * nNH3 * nET * xC**2 * xCH2**2 * xET**2 * xHe**2 * a**2 * \
    #         (1.0 * nHe * xC * xCH2 - 1.19655 * nCH2 * xC * xHe + 0.196547 * nC * xCH2 * xHe)**2
    # term2 = 85044.9 * nNH3 * nCH2 * xC**2 * xCH2**2 * xET**2 * xHe**2 * a**2 * \
    #         (1.0 * nHe * xC * xET - 0.302592 * nET * xC * xHe - 0.697408 * nC * xET * xHe)**2
    # term3 = 47470.0 * nNH3 * nC * xC**2 * xCH2**2 * xET**2 * xHe**2 * a**2 * \
    #         (1.0 * nHe * xCH2 * xET - 0.0665285 * nET * xCH2 * xHe - 0.933472 * nCH2 * xET * xHe)**2
    # term4 = 1371.83 * nHe**2 * xNH3**2 * \
    #         (1.0 * nHe * xC * xCH2 * xET + (-1.97748 * nET * xC * xCH2 + 6.62306 * nCH2 * xC * xET - 5.64558 * nC * xCH2 * xET) * xHe)**2 * \
    #         (1.0 * nHe * xC * xCH2 * xET + (0.0136533 * nET * xC * xCH2 - 1.25054 * nCH2 * xC * xET + 0.236883 * nC * xCH2 * xET) * xHe)**2
    # t5_p1 = 73.2426 * nHe**2 * xNH3 * xC**2 * xCH2**2 * xET**2
    # t5_p2 = nHe * xC * xCH2 * xET * \
    #         (2.0 * nET * xNH3 * xC * xCH2 + (-183.185 * nCH2 * xNH3 * xC + 34.6998 * nC * xNH3 * xCH2 + 1.42109e-14 * nNH3 * xC * xCH2) * xET)
    # t5_p3 = (-1.97748 * nET**2 * xNH3 * xC**2 * xCH2**2) + \
    #          nET * xC * xCH2 * (187.746 * nCH2 * xNH3 * xC - 39.9547 * nC * xNH3 * xCH2 - 145.836 * nNH3 * xC * xCH2) * xET + \
    #          (-606.623 * nCH2**2 * xNH3 * xC**2 + nCH2 * xC * (632.002 * nC * xNH3 + 576.683 * nNH3 * xC) * xCH2 + \
    #          nC * (-97.9502 * nC * xNH3 - 430.847 * nNH3 * xC) * xCH2**2) * xET**2
    # term5 = 0.255725 * nNH3 * nHe * (t5_p1 + t5_p2 * xHe + t5_p3 * xHe**2)**2
    # denom = nNH3**3 * xHe**2 * \
    #         (73.2426 * nHe * xC * xCH2 * xET + 1.0 * nET * xC * xCH2 * xHe - 91.5925 * nCH2 * xC * xET * xHe + 17.3499 * nC * xCH2 * xET * xHe)**4
    # total = term1 + term2 + term3 + term4 + term5
    # print(f"[Debug Error] sum_terms={total:.4e}, denom={denom:.4e}")
    # var = total/denom if denom else 0.0
    # if var < 0:
    #     print(f"[Warning] Negative var={var:.4e}, clamping to 0.")
    #     var = 0.0
    # sigma_df = 27.3473 * math.sqrt(var)
    # First part of the expression (term1)
    double term1 = 5438.81 * nNH3 * nET * pow(xC, 2) * pow(xCH2, 2) * pow(xET, 2) * pow(xHe, 2)
                   * pow((-1.0 * nHe * xNH3 + nNH3 * xHe), 2)
                   * pow((1.0 * nHe * xC * xCH2 - 1.19655 * nCH2 * xC * xHe + 0.196547 * nC * xCH2 * xHe), 2);

    # Second part of the expression (term2)
    double term2 = 85044.9 * nNH3 * nCH2 * pow(xC, 2) * pow(xCH2, 2) * pow(xET, 2) * pow(xHe, 2)
                   * pow((-1.0 * nHe * xNH3 + nNH3 * xHe), 2)
                   * pow((1.0 * nHe * xC * xET - 0.302592 * nET * xC * xHe - 0.697408 * nC * xET * xHe), 2);

    # Third part of the expression (term3)
    double term3 = 47470.0 * nNH3 * nC * pow(xC, 2) * pow(xCH2, 2) * pow(xET, 2) * pow(xHe, 2)
                   * pow((-1.0 * nHe * xNH3 + nNH3 * xHe), 2)
                   * pow((1.0 * nHe * xCH2 * xET - 0.0665285 * nET * xCH2 * xHe - 0.933472 * nCH2 * xET * xHe), 2);

    # Fourth part of the expression (term4)
    double term4 = 1371.83 * pow(nHe, 2) * pow(xNH3, 2)
                   * pow((1.0 * nHe * xC * xCH2 * xET + (-1.97748 * nET * xC * xCH2 + 6.62306 * nCH2 * xC * xET - 5.64558 * nC * xCH2 * xET) * xHe), 2)
                   * pow((1.0 * nHe * xC * xCH2 * xET + (0.0136533 * nET * xC * xCH2 - 1.25054 * nCH2 * xC * xET + 0.236883 * nC * xCH2 * xET) * xHe), 2);

    # Fifth part of the expression (term5)
    # First part of the numerator
    double term5_numerator_part1 = 73.2426 * pow(nHe, 2) * xNH3 * pow(xC, 2) * pow(xCH2, 2) * pow(xET, 2);

    # Second part of the numerator (to be multiplied by xHe)
    double term5_numerator_part2 = nHe * xC * xCH2 * xET
                                   * (2.0 * nET * xNH3 * xC * xCH2 + (-183.185 * nCH2 * xNH3 * xC + 34.6998 * nC * xNH3 * xCH2 + 1.42109e-14 * nNH3 * xC * xCH2) * xET);

    # Third part of the numerator (to be multiplied by xHe^2)
    double term5_numerator_part3 = (-1.97748 * pow(nET, 2) * xNH3 * pow(xC, 2) * pow(xCH2, 2))
                                   + nET * xC * xCH2 * (187.746 * nCH2 * xNH3 * xC - 39.9547 * nC * xNH3 * xCH2 - 145.836 * nNH3 * xC * xCH2) * xET
                                   + (-606.623 * pow(nCH2, 2) * xNH3 * pow(xC, 2)
                                      + nCH2 * xC * (632.002 * nC * xNH3 + 576.683 * nNH3 * xC) * xCH2
                                      + nC * (-97.9502 * nC * xNH3 - 430.847 * nNH3 * xC) * pow(xCH2, 2)
                                     ) * pow(xET, 2);

    # The full numerator
    double term5_numerator = term5_numerator_part1 + term5_numerator_part2 * xHe + term5_numerator_part3 * pow(xHe, 2);

    double term5 = 0.255725 * nNH3 * nHe * pow(term5_numerator, 2);

    # DenominNH3tor of the expression
    double denominNH3tor = pow(nNH3, 3) * pow(xHe, 2)
                         * pow((73.2426 * nHe * xC * xCH2 * xET + 1.0 * nET * xC * xCH2 * xHe - 91.5925 * nCH2 * xC * xET * xHe + 17.3499 * nC * xCH2 * xET * xHe), 4);

    # FinNH3l error calculation
    double sigma_df = 27.3473 * sqrt((term1 + term2 + term3 + term4 + term5) / denominNH3tor);

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
        # mask y
        idx = {t: np.digitize(data_x[t][data_y[t]<Y_CUT], xB_bins)-1 for t in data_x}
        x_mean = np.mean(data_x["NH3"][data_y["NH3"]<Y_CUT][idx["NH3"]==i]) if np.any(idx["NH3"]==i) else 0.0
        counts = {t: np.sum(idx[t]==i) for t in idx}
        print(f"[Period {period}] Bin{i+1}/{nbins} x_mean={x_mean:.4f} counts={counts}")
        Df, Df_s = calculate_dilution_and_error(
            counts["NH3"],counts["C"],counts["CH2"],counts["He"],counts["ET"],
            fr["xNH3"],fr["xC"],fr["xCH2"],fr["xHe"],fr["xET"]
        )
        rows.append({"x_mean":x_mean,"Df":Df,"Df_sigma":Df_s})
    df = pd.DataFrame(rows); df["period"]=period
    print(f"[Period {period}] Done.\n")
    return df


def calculate_and_save(trees, xB_bins):
    print("[Main] Begin dilution calculation")
    os.makedirs("output",exist_ok=True)
    per = ["RGC_Su22","RGC_Fa22"]
    with ProcessPoolExecutor() as ex:
        futs={ex.submit(_compute_period,p,trees[p],xB_bins):p for p in per}
        res={futs[f]:f.result() for f in as_completed(futs)}
    df_su,df_fa=res["RGC_Su22"],res["RGC_Fa22"]
    wsu=1/df_su["Df_sigma"]**2; wfa=1/df_fa["Df_sigma"]**2
    df_avg=pd.DataFrame({"x_mean":df_su["x_mean"],"Df":(df_su["Df"]*wsu+df_fa["Df"]*wfa)/(wsu+wfa),"Df_sigma":np.sqrt(1/(wsu+wfa)),"period":"weighted"})
    df_all=pd.concat([df_su,df_fa,df_avg],ignore_index=True)
    csv="output/dilution_factor.csv"; df_all.to_csv(csv,index=False); print(f"Saved CSV {csv}")
    plt.errorbar(df_su["x_mean"],df_su["Df"],yerr=df_su["Df_sigma"],fmt="o",label="Su22")
    plt.errorbar(df_fa["x_mean"],df_fa["Df"],yerr=df_fa["Df_sigma"],fmt="o",label="Fa22")
    plt.xlabel("$x_B$");plt.ylabel("$D_f$");plt.legend(loc="upper right");plt.tight_layout()
    pdf="output/dilution_factor.pdf"; plt.savefig(pdf); plt.close(); print(f"Saved PDF {pdf}\nDone.")
