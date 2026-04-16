#!/usr/bin/env python3

import os
import csv
import math
import uuid
import numpy as np
import matplotlib.pyplot as plt
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# ------------------------------------------------
# configuration
# ------------------------------------------------

INPUT_FILE = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_combined_inb_NH3_epi+_2.root"
TREE_NAME = "PhysicsEvents"
OUTPUT_DIR = "output/enpi+_Mx2_fits"

# histogram / fit settings
MX2_MIN = 0.70
MX2_MAX = 1.10
N_MX2_BINS = 60

MU_MIN = 0.85
MU_MAX = 0.95
SIGMA_MIN = 0.03
SIGMA_MAX = 0.08

# same numerical binning scheme as before, but applied to positive (-tprime)
XB_SLICES = [
    ("Low", 0.10, 0.25),
    ("MidLow", 0.25, 0.35),
    ("MidHigh", 0.35, 0.45),
    ("High", 0.45, 0.60),
]

TPRIME_BINS_POS = [
    (0.05, 0.25),
    (0.25, 0.45),
    (0.45, 0.65),
    (0.65, 0.85),
    (0.85, 1.05),
    (1.05, 1.25),
]

CSV_PATH = os.path.join(OUTPUT_DIR, "Mx2_fit_params_tprime.csv")
LATEX_TABLE_PATH = os.path.join(OUTPUT_DIR, "Mx2_fit_params_tprime_table.tex")
GRID_PNG_PATH = os.path.join(OUTPUT_DIR, "Mx2_tprime_fit_grid.png")
GRID_PDF_PATH = os.path.join(OUTPUT_DIR, "Mx2_tprime_fit_grid.pdf")
SUMMARY_PNG_PATH = os.path.join(OUTPUT_DIR, "Mx2_tprime_mu_sigma_summary.png")
SUMMARY_PDF_PATH = os.path.join(OUTPUT_DIR, "Mx2_tprime_mu_sigma_summary.pdf")

# ------------------------------------------------
# helpers
# ------------------------------------------------

def ensure_output_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)
    #endif
#enddef

def load_arrays(root_file, tree_name):
    rdf = ROOT.RDataFrame(tree_name, root_file)
    arrays = rdf.AsNumpy(["Q2", "W", "y", "x", "tprime", "Mx2"])
    out = {}
    for key in arrays:
        out[key] = np.asarray(arrays[key], dtype=np.float64)
    #endfor
    return out
#enddef

def base_event_mask(arrays):
    mask = (
        (arrays["Q2"] > 1.0) &
        (arrays["W"] > 2.0) &
        (arrays["y"] < 0.75) &
        np.isfinite(arrays["x"]) &
        np.isfinite(arrays["tprime"]) &
        np.isfinite(arrays["Mx2"])
    )
    return mask
#enddef

def make_histogram(values, xmin, xmax, nbins):
    counts, edges = np.histogram(values, bins=nbins, range=(xmin, xmax))
    centers = 0.5 * (edges[:-1] + edges[1:])
    errors = np.sqrt(np.maximum(counts.astype(np.float64), 1.0))
    return counts.astype(np.float64), edges, centers, errors
#enddef

def build_root_hist(counts, edges, name):
    hist = ROOT.TH1D(name, "", len(counts), float(edges[0]), float(edges[-1]))
    for i in range(len(counts)):
        hist.SetBinContent(i + 1, float(counts[i]))
        hist.SetBinError(i + 1, math.sqrt(max(float(counts[i]), 1.0)))
    #endfor
    return hist
#enddef

def clamp(value, vmin, vmax):
    return max(vmin, min(vmax, value))
#enddef

def fit_mx2_hist(counts, edges):
    hist_name = "h_" + str(uuid.uuid4()).replace("-", "_")
    func_name = "f_" + str(uuid.uuid4()).replace("-", "_")

    hist = build_root_hist(counts, edges, hist_name)

    peak_bin = int(np.argmax(counts))
    peak_x = 0.5 * (edges[peak_bin] + edges[peak_bin + 1])
    edge_background = 0.5 * (np.mean(counts[:5]) + np.mean(counts[-5:]))
    amplitude_guess = max(float(np.max(counts) - edge_background), 1.0)
    mu_guess = clamp(float(peak_x), MU_MIN + 0.002, MU_MAX - 0.002)
    sigma_guess = 0.07

    fit_fn = ROOT.TF1(func_name, "gaus(0) + pol2(3)", MX2_MIN, MX2_MAX)
    fit_fn.SetParameter(0, amplitude_guess)
    fit_fn.SetParameter(1, mu_guess)
    fit_fn.SetParameter(2, sigma_guess)
    fit_fn.SetParameter(3, edge_background)
    fit_fn.SetParameter(4, 0.0)
    fit_fn.SetParameter(5, 0.0)

    fit_fn.SetParLimits(0, 0.0, 1.0e12)
    fit_fn.SetParLimits(1, MU_MIN, MU_MAX)
    fit_fn.SetParLimits(2, SIGMA_MIN, SIGMA_MAX)

    fit_result = hist.Fit(fit_fn, "QSRN")

    status = int(fit_result.Status())
    cov_status = int(fit_result.CovMatrixStatus())
    fit_success = int(status == 0 and cov_status >= 2)

    result = {
        "fit_success": fit_success,
        "status": status,
        "cov_status": cov_status,
        "chi2": float(fit_fn.GetChisquare()),
        "ndf": int(fit_fn.GetNDF()),
        "amp": float("nan"),
        "mu": float("nan"),
        "mu_err": float("nan"),
        "sigma": float("nan"),
        "sigma_err": float("nan"),
        "b0": float("nan"),
        "b1": float("nan"),
        "b2": float("nan"),
    }

    if fit_success:
        result["amp"] = float(fit_fn.GetParameter(0))
        result["mu"] = float(fit_fn.GetParameter(1))
        result["mu_err"] = float(fit_fn.GetParError(1))
        result["sigma"] = float(fit_fn.GetParameter(2))
        result["sigma_err"] = float(fit_fn.GetParError(2))
        result["b0"] = float(fit_fn.GetParameter(3))
        result["b1"] = float(fit_fn.GetParameter(4))
        result["b2"] = float(fit_fn.GetParameter(5))
    #endif

    ROOT.SetOwnership(hist, False)
    ROOT.SetOwnership(fit_fn, False)

    return result
#enddef

def evaluate_fit_curve(xvals, fit_result):
    amp = fit_result["amp"]
    mu = fit_result["mu"]
    sigma = fit_result["sigma"]
    b0 = fit_result["b0"]
    b1 = fit_result["b1"]
    b2 = fit_result["b2"]

    gaus = amp * np.exp(-0.5 * ((xvals - mu) / sigma) ** 2)
    poly = b0 + b1 * xvals + b2 * xvals * xvals
    return gaus + poly
#enddef

def make_grid_figure(all_rows):
    fig, axes = plt.subplots(
        nrows=len(XB_SLICES),
        ncols=len(TPRIME_BINS_POS),
        figsize=(18, 12),
        sharex=True,
        constrained_layout=False
    )

    fig.patch.set_facecolor("#efefef")

    for row in range(len(XB_SLICES)):
        for col in range(len(TPRIME_BINS_POS)):
            ax = axes[row, col]
            ax.set_facecolor("#efefef")
            ax.grid(True, linestyle="--", alpha=0.35)
            ax.tick_params(axis="both", labelsize=8)

            entry = all_rows[row * len(TPRIME_BINS_POS) + col]

            centers = entry["centers"]
            counts = entry["counts"]
            errors = entry["errors"]

            ax.errorbar(
                centers,
                counts,
                yerr=errors,
                fmt="s",
                markersize=2.5,
                linewidth=0.8,
                elinewidth=0.8,
                capsize=1.2,
                color="black"
            )

            if entry["fit_success"] == 1:
                xfit = np.linspace(MX2_MIN, MX2_MAX, 600)
                yfit = evaluate_fit_curve(xfit, entry)
                ax.plot(
                    xfit,
                    yfit,
                    linestyle="--",
                    linewidth=1.2,
                    color="red",
                    label=r"$\mu = %.4f,\ \sigma = %.4f$" % (entry["mu"], entry["sigma"])
                )
            else:
                ax.plot([], [], linestyle="--", linewidth=1.2, color="red", label="fit failed")
            #endif

            ax.legend(loc="lower center", fontsize=7, framealpha=0.9)

            xlo = entry["x_min"]
            xhi = entry["x_max"]
            tlo = entry["t_min"]
            thi = entry["t_max"]

            ax.set_title(
                r"$x_B \in [%.2f, %.2f],\ -t^\prime \in [%.2f, %.2f]$" % (xlo, xhi, tlo, thi),
                fontsize=8
            )

            ax.set_xlim(MX2_MIN, MX2_MAX)

            if col == 0:
                ax.set_ylabel("counts", fontsize=9)
            #endif

            if row == len(XB_SLICES) - 1:
                ax.set_xlabel(r"$M_X^2$ (GeV$^2$)", fontsize=9)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        r"$M_X^2$ distributions by $x_B$ slice and $-t^\prime$ bin" + "\n" +
        r"Global cuts: $Q^2 > 1$, $W > 2$, $y < 0.75$ (no $M_X^2$ cut)",
        fontsize=16,
        y=0.98
    )

    plt.subplots_adjust(left=0.06, right=0.995, bottom=0.07, top=0.90, wspace=0.12, hspace=0.15)
    fig.savefig(GRID_PNG_PATH, dpi=300)
    fig.savefig(GRID_PDF_PATH)
    plt.close(fig)
#enddef

def write_csv(rows, csv_path):
    fieldnames = [
        "slice_name",
        "x_min",
        "x_max",
        "t_min",
        "t_max",
        "mu",
        "mu_err",
        "sigma",
        "sigma_err",
        "n_entries",
        "fit_success",
    ]

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "slice_name": row["slice_name"],
                "x_min": f"{row['x_min']:.4f}",
                "x_max": f"{row['x_max']:.4f}",
                "t_min": f"{row['t_min']:.4f}",
                "t_max": f"{row['t_max']:.4f}",
                "mu": "" if not np.isfinite(row["mu"]) else f"{row['mu']:.6f}",
                "mu_err": "" if not np.isfinite(row["mu_err"]) else f"{row['mu_err']:.6f}",
                "sigma": "" if not np.isfinite(row["sigma"]) else f"{row['sigma']:.6f}",
                "sigma_err": "" if not np.isfinite(row["sigma_err"]) else f"{row['sigma_err']:.6f}",
                "n_entries": int(row["n_entries"]),
                "fit_success": int(row["fit_success"]),
            })
        #endfor
    #endwith
#enddef

def write_latex_table(rows, latex_path):
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering")
    lines.append(r"  \caption{Fit results in $(x_B,-t^\prime)$ bins. Shown are the fitted mean $\mu$, its uncertainty $\delta\mu$, the fitted width $\sigma$, and its uncertainty $\delta\sigma$.}")
    lines.append(r"  \label{tab:mx2_tprime_fit_results}")
    lines.append(r"  \begin{tabular}{cccccccc}")
    lines.append(r"    \hline")
    lines.append(r"    $x_{B,\min}$ & $x_{B,\max}$ & $(-t^\prime)_{\min}$ & $(-t^\prime)_{\max}$ & $\mu$ & $\delta\mu$ & $\sigma$ & $\delta\sigma$ \\")
    lines.append(r"    \hline")

    last_slice = None
    for row in rows:
        if last_slice is not None and row["slice_name"] != last_slice:
            lines.append(r"    \hline")
        #endif

        mu_str = "--"
        mu_err_str = "--"
        sigma_str = "--"
        sigma_err_str = "--"

        if row["fit_success"] == 1 and np.isfinite(row["mu"]):
            mu_str = f"{row['mu']:.4f}"
            mu_err_str = f"{row['mu_err']:.4f}"
            sigma_str = f"{row['sigma']:.4f}"
            sigma_err_str = f"{row['sigma_err']:.4f}"
        #endif

        lines.append(
            "    "
            + f"{row['x_min']:.4f} & "
            + f"{row['x_max']:.4f} & "
            + f"{row['t_min']:.4f} & "
            + f"{row['t_max']:.4f} & "
            + f"{mu_str} & "
            + f"{mu_err_str} & "
            + f"{sigma_str} & "
            + f"{sigma_err_str} \\\\"
        )

        last_slice = row["slice_name"]
    #endfor

    lines.append(r"    \hline")
    lines.append(r"  \end{tabular}")
    lines.append(r"\end{table}")

    with open(latex_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    #endwith
#enddef

def make_summary_figure(rows):
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(14, 9),
        sharex=True,
        constrained_layout=False
    )

    fig.patch.set_facecolor("#efefef")
    colors = ["C0", "C1", "C2", "C3"]
    labels = [
        r"$x_B \in [0.10, 0.25]$",
        r"$x_B \in [0.25, 0.35]$",
        r"$x_B \in [0.35, 0.45]$",
        r"$x_B \in [0.45, 0.60]$",
    ]

    x_positions = np.arange(1, len(rows) + 1)

    for ax in axes:
        ax.set_facecolor("#efefef")
        ax.grid(True, linestyle="--", alpha=0.45)
        ax.tick_params(axis="both", labelsize=11)
        ax.axvline(6.5, color="0.75", linestyle=":", linewidth=1.0)
        ax.axvline(12.5, color="0.75", linestyle=":", linewidth=1.0)
        ax.axvline(18.5, color="0.75", linestyle=":", linewidth=1.0)
    #endfor

    for i_slice in range(len(XB_SLICES)):
        start = i_slice * len(TPRIME_BINS_POS)
        stop = start + len(TPRIME_BINS_POS)
        subrows = rows[start:stop]

        xpos = x_positions[start:stop]
        mu = np.array([r["mu"] for r in subrows], dtype=np.float64)
        mu_err = np.array([r["mu_err"] for r in subrows], dtype=np.float64)
        sigma = np.array([r["sigma"] for r in subrows], dtype=np.float64)
        sigma_err = np.array([r["sigma_err"] for r in subrows], dtype=np.float64)

        axes[0].errorbar(
            xpos,
            mu,
            yerr=mu_err,
            fmt="o",
            markersize=5,
            linewidth=1.5,
            capsize=4,
            color=colors[i_slice],
            label=labels[i_slice]
        )

        axes[1].errorbar(
            xpos,
            sigma,
            yerr=sigma_err,
            fmt="o",
            markersize=5,
            linewidth=1.5,
            capsize=4,
            color=colors[i_slice],
            label=labels[i_slice]
        )
    #endfor

    axes[0].set_ylabel(r"$\mu$", fontsize=14)
    axes[1].set_ylabel(r"$\sigma$", fontsize=14)
    axes[1].set_xlabel("bin", fontsize=14)

    finite_mu = np.array([r["mu"] for r in rows if np.isfinite(r["mu"])], dtype=np.float64)
    finite_mu_err = np.array([r["mu_err"] for r in rows if np.isfinite(r["mu_err"])], dtype=np.float64)
    finite_sigma = np.array([r["sigma"] for r in rows if np.isfinite(r["sigma"])], dtype=np.float64)
    finite_sigma_err = np.array([r["sigma_err"] for r in rows if np.isfinite(r["sigma_err"])], dtype=np.float64)

    if len(finite_mu) > 0:
        mu_lo = np.min(finite_mu - finite_mu_err) - 0.01
        mu_hi = np.max(finite_mu + finite_mu_err) + 0.01
        axes[0].set_ylim(mu_lo, mu_hi)
    #endif

    if len(finite_sigma) > 0:
        sigma_lo = max(0.0, np.min(finite_sigma - finite_sigma_err) - 0.005)
        sigma_hi = np.max(finite_sigma + finite_sigma_err) + 0.005
        axes[1].set_ylim(sigma_lo, sigma_hi)
    #endif

    axes[0].legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.18),
        ncol=4,
        fontsize=12,
        framealpha=1.0
    )

    plt.subplots_adjust(left=0.08, right=0.99, bottom=0.08, top=0.88, hspace=0.14)
    fig.savefig(SUMMARY_PNG_PATH, dpi=300)
    fig.savefig(SUMMARY_PDF_PATH)
    plt.close(fig)
#enddef

# ------------------------------------------------
# main analysis
# ------------------------------------------------

def main():
    ensure_output_dir(OUTPUT_DIR)

    arrays = load_arrays(INPUT_FILE, TREE_NAME)
    mask_base = base_event_mask(arrays)

    x_all = arrays["x"][mask_base]
    tprime_all = arrays["tprime"][mask_base]
    mx2_all = arrays["Mx2"][mask_base]

    # tprime in the tree is negative, so the binned quantity is -tprime
    tprime_pos_all = -tprime_all

    all_rows = []

    for slice_name, x_min, x_max in XB_SLICES:
        for t_min, t_max in TPRIME_BINS_POS:
            sel = (
                (x_all > x_min) &
                (x_all < x_max) &
                (tprime_pos_all >= t_min) &
                (tprime_pos_all < t_max)
            )

            mx2_sel = mx2_all[sel]
            counts, edges, centers, errors = make_histogram(mx2_sel, MX2_MIN, MX2_MAX, N_MX2_BINS)
            fit = fit_mx2_hist(counts, edges)

            row = {
                "slice_name": slice_name,
                "x_min": x_min,
                "x_max": x_max,
                "t_min": t_min,
                "t_max": t_max,
                "n_entries": int(len(mx2_sel)),
                "fit_success": int(fit["fit_success"]),
                "status": int(fit["status"]),
                "cov_status": int(fit["cov_status"]),
                "chi2": float(fit["chi2"]),
                "ndf": int(fit["ndf"]),
                "mu": float(fit["mu"]),
                "mu_err": float(fit["mu_err"]),
                "sigma": float(fit["sigma"]),
                "sigma_err": float(fit["sigma_err"]),
                "amp": float(fit["amp"]),
                "b0": float(fit["b0"]),
                "b1": float(fit["b1"]),
                "b2": float(fit["b2"]),
                "counts": counts,
                "edges": edges,
                "centers": centers,
                "errors": errors,
            }

            all_rows.append(row)

            print(
                f"{slice_name:8s}  "
                f"x=[{x_min:.2f},{x_max:.2f}]  "
                f"-tprime=[{t_min:.2f},{t_max:.2f}]  "
                f"N={len(mx2_sel):7d}  "
                f"fit_success={row['fit_success']}  "
                f"mu={row['mu']:.6f}  "
                f"mu_err={row['mu_err']:.6f}  "
                f"sigma={row['sigma']:.6f}  "
                f"sigma_err={row['sigma_err']:.6f}"
            )
        #endfor
    #endfor

    write_csv(all_rows, CSV_PATH)
    write_latex_table(all_rows, LATEX_TABLE_PATH)
    make_grid_figure(all_rows)
    make_summary_figure(all_rows)

    print("")
    print("Done.")
    print(f"CSV:          {CSV_PATH}")
    print(f"LaTeX table:  {LATEX_TABLE_PATH}")
    print(f"Grid figure:  {GRID_PNG_PATH}")
    print(f"Summary fig:  {SUMMARY_PNG_PATH}")
#enddef

if __name__ == "__main__":
    main()
#endif