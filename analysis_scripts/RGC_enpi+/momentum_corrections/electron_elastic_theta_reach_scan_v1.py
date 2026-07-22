#!/usr/bin/env python3
"""
Scan how far in electron polar angle the inclusive elastic eX sample can
support a usable W-peak calibration.

For each RGC period and FD sector, this script:
  * fills fixed-width electron-theta bins;
  * fits W with Gaussian + quadratic background;
  * records event count, fitted centroid, centroid uncertainty, sigma,
    peak significance, and chi2/ndf;
  * classifies each theta bin as usable or unusable;
  * reports the highest contiguous usable theta bin;
  * makes coverage and fit-quality plots.

This is a diagnostic scan only. It does not change the production correction.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from pathlib import Path
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit


PROTON_MASS_GEV = 0.9382720813
RAD2DEG = 180.0 / math.pi

DEFAULT_INPUTS = {
    "su22": Path(
        "/work/clas12/thayward/CLAS12_exclusive/elastic/"
        "rgc_su22_inb_eX_Wlessthan2.root"
    ),
    "fa22": Path(
        "/work/clas12/thayward/CLAS12_exclusive/elastic/"
        "rgc_fa22_inb_eX_Wlessthan2.root"
    ),
    "sp23": Path(
        "/work/clas12/thayward/CLAS12_exclusive/elastic/"
        "rgc_sp23_inb_eX_Wlessthan2.root"
    ),
}


@dataclass(frozen=True)
class Period:
    key: str
    label: str
    source_key: str
    run_min: int
    run_max: int


PERIODS = (
    Period("su22", "Su22", "su22", 16042, 16788),
    Period("fa22_solenoid_minus", "Fa22, solenoid -1", "fa22", 16843, 17183),
    Period("fa22_solenoid_plus", "Fa22, solenoid +1", "fa22", 17185, 17408),
    Period("sp23", "Sp23", "sp23", 17477, 17811),
)


def fd_sector_from_phi(phi_rad: np.ndarray) -> np.ndarray:
    phi_deg = np.mod(np.degrees(phi_rad), 360.0)
    sector = np.full(phi_deg.shape, -1, dtype=np.int16)
    sector[(phi_deg >= 330.0) | (phi_deg < 30.0)] = 1
    sector[(phi_deg >= 30.0) & (phi_deg < 90.0)] = 2
    sector[(phi_deg >= 90.0) & (phi_deg < 150.0)] = 3
    sector[(phi_deg >= 150.0) & (phi_deg < 210.0)] = 4
    sector[(phi_deg >= 210.0) & (phi_deg < 270.0)] = 5
    sector[(phi_deg >= 270.0) & (phi_deg < 330.0)] = 6
    return sector


def model(x, amplitude, mean, sigma, b0, b1, b2):
    dx = x - PROTON_MASS_GEV
    return (
        amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
        + b0 + b1 * dx + b2 * dx**2
    )


def fit_one_bin(
    w: np.ndarray,
    w_hist_min: float,
    w_hist_max: float,
    w_fit_min: float,
    w_fit_max: float,
    w_bins: int,
) -> dict:
    counts, edges = np.histogram(
        w, bins=w_bins, range=(w_hist_min, w_hist_max)
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    mask = (centers >= w_fit_min) & (centers <= w_fit_max)
    x = centers[mask]
    y = counts[mask].astype(float)
    ey = np.sqrt(np.maximum(y, 1.0))

    if np.count_nonzero(y) < 10:
        return {"fit_success": False, "fit_status": "too_few_populated_bins"}

    peak_mask = (x >= 0.88) & (x <= 1.00)
    if not np.any(peak_mask):
        return {"fit_success": False, "fit_status": "empty_peak_window"}

    sideband = (x < 0.88) | (x > 1.00)
    if np.count_nonzero(sideband) >= 6:
        try:
            coeff = np.polyfit(
                x[sideband] - PROTON_MASS_GEV,
                y[sideband],
                2,
                w=1.0 / ey[sideband],
            )
            bg_seed = np.polyval(coeff, x - PROTON_MASS_GEV)
            b2, b1, b0 = map(float, coeff)
        except Exception:
            bg_seed = np.full_like(y, np.median(y[sideband]))
            b0, b1, b2 = float(np.median(y[sideband])), 0.0, 0.0
    else:
        bg_seed = np.full_like(y, np.percentile(y, 25))
        b0, b1, b2 = float(np.percentile(y, 25)), 0.0, 0.0

    residual = y - bg_seed
    local_indices = np.flatnonzero(peak_mask)
    peak_index = local_indices[np.argmax(residual[peak_mask])]

    p0 = (
        max(float(residual[peak_index]), 1.0),
        float(np.clip(x[peak_index], 0.88, 1.00)),
        0.025,
        b0, b1, b2,
    )
    lower = (0.0, 0.86, 0.008, -np.inf, -np.inf, -np.inf)
    upper = (np.inf, 1.02, 0.080, np.inf, np.inf, np.inf)

    try:
        pars, cov = curve_fit(
            model, x, y, p0=p0, sigma=ey, absolute_sigma=True,
            bounds=(lower, upper), maxfev=100000
        )
    except Exception as exc:
        return {
            "fit_success": False,
            "fit_status": f"fit_failed:{type(exc).__name__}",
        }

    errors = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    amplitude, mean, sigma, _, _, _ = pars
    amplitude_error, mean_error, sigma_error = errors[:3]
    prediction = model(x, *pars)
    chi2 = float(np.sum(((y - prediction) / ey) ** 2))
    ndf = int(len(x) - len(pars))
    chi2_ndf = chi2 / ndf if ndf > 0 else math.nan
    significance = (
        float(amplitude / amplitude_error)
        if np.isfinite(amplitude_error) and amplitude_error > 0
        else math.nan
    )

    return {
        "fit_success": True,
        "fit_status": "success",
        "peak_mean_gev": float(mean),
        "peak_mean_error_gev": float(mean_error),
        "peak_sigma_gev": float(abs(sigma)),
        "peak_sigma_error_gev": float(sigma_error),
        "amplitude": float(amplitude),
        "amplitude_error": float(amplitude_error),
        "peak_significance": significance,
        "chi2": chi2,
        "ndf": ndf,
        "chi2_ndf": chi2_ndf,
    }


def scan_period(period: Period, inputs: dict[str, Path], args) -> pd.DataFrame:
    path = inputs[period.source_key]
    source = f"{path}:{args.tree_name}"

    chunks = []
    for arrays in uproot.iterate(
        source,
        expressions=["runnum", "e_theta", "e_phi", "W"],
        library="np",
        step_size=args.step_size,
    ):
        run = np.asarray(arrays["runnum"], dtype=np.int64)
        theta = np.degrees(np.asarray(arrays["e_theta"], dtype=float))
        phi = np.asarray(arrays["e_phi"], dtype=float)
        w = np.asarray(arrays["W"], dtype=float)

        mask = (
            (run >= period.run_min) & (run <= period.run_max)
            & np.isfinite(theta) & np.isfinite(phi) & np.isfinite(w)
            & (theta >= args.theta_min) & (theta < args.theta_max)
            & (w >= args.w_preselection_min)
            & (w <= args.w_preselection_max)
        )
        if not np.any(mask):
            continue

        sector = fd_sector_from_phi(phi[mask])
        good = (sector >= 1) & (sector <= 6)
        if not np.any(good):
            continue

        chunks.append(pd.DataFrame({
            "theta_deg": theta[mask][good],
            "sector": sector[good],
            "W": w[mask][good],
        }))

    if not chunks:
        raise RuntimeError(f"{period.label}: no selected events")

    data = pd.concat(chunks, ignore_index=True)
    edges = np.arange(
        args.theta_min,
        args.theta_max + 0.5 * args.theta_bin_width,
        args.theta_bin_width,
    )

    records = []
    for sector in range(1, 7):
        sector_data = data[data["sector"] == sector]
        theta_values = sector_data["theta_deg"].to_numpy()
        w_values = sector_data["W"].to_numpy()

        for low, high in zip(edges[:-1], edges[1:]):
            mask = (theta_values >= low) & (theta_values < high)
            selected_w = w_values[mask]
            record = {
                "period_key": period.key,
                "period": period.label,
                "sector": sector,
                "theta_low_deg": float(low),
                "theta_high_deg": float(high),
                "theta_center_deg": float(0.5 * (low + high)),
                "n_events": int(selected_w.size),
            }

            if selected_w.size < args.minimum_events:
                record.update({
                    "fit_success": False,
                    "fit_status": "insufficient_events",
                })
            else:
                record.update(fit_one_bin(
                    selected_w,
                    args.w_hist_min,
                    args.w_hist_max,
                    args.w_fit_min,
                    args.w_fit_max,
                    args.w_bins,
                ))

            records.append(record)

    return pd.DataFrame(records)


def classify(df: pd.DataFrame, args) -> pd.DataFrame:
    result = df.copy()

    for column in (
        "peak_mean_gev", "peak_mean_error_gev", "peak_sigma_gev",
        "peak_significance", "chi2_ndf"
    ):
        if column not in result:
            result[column] = np.nan

    result["usable"] = (
        result["fit_success"].fillna(False)
        & (result["n_events"] >= args.minimum_events)
        & (result["peak_mean_error_gev"] <= args.max_mean_error)
        & (result["peak_significance"] >= args.minimum_significance)
        & (result["chi2_ndf"] <= args.max_chi2_ndf)
        & (result["peak_sigma_gev"] >= args.min_sigma)
        & (result["peak_sigma_gev"] <= args.max_sigma)
        & (result["peak_mean_gev"] >= args.min_peak_mean)
        & (result["peak_mean_gev"] <= args.max_peak_mean)
    )
    return result


def contiguous_reach(group: pd.DataFrame) -> float:
    ordered = group.sort_values("theta_low_deg")
    usable = ordered["usable"].to_numpy(dtype=bool)
    if not np.any(usable):
        return math.nan

    # Begin at the first usable bin, then stop at the first subsequent failure.
    start = int(np.argmax(usable))
    last = start
    for index in range(start + 1, len(usable)):
        if not usable[index]:
            break
        last = index
    return float(ordered.iloc[last]["theta_high_deg"])


def make_plots(df: pd.DataFrame, output_dir: Path, args) -> None:
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    for period, period_df in df.groupby("period", sort=False):
        safe = (
            period.lower()
            .replace(",", "")
            .replace(" ", "_")
            .replace("+", "plus")
            .replace("-", "minus")
        )

        fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True)
        axes = axes.ravel()
        for sector in range(1, 7):
            ax = axes[sector - 1]
            sdf = period_df[period_df["sector"] == sector]
            ax.step(
                sdf["theta_center_deg"],
                sdf["n_events"],
                where="mid",
                label="Selected entries",
            )
            bad = sdf[~sdf["usable"]]
            good = sdf[sdf["usable"]]
            ax.scatter(good["theta_center_deg"], good["n_events"], marker="o",
                       label="Usable fit")
            ax.scatter(bad["theta_center_deg"], bad["n_events"], marker="x",
                       label="Rejected / unavailable")
            ax.set_yscale("log")
            ax.set_title(f"Sector {sector}")
            ax.set_xlabel(r"$\theta_e$ (deg)")
            ax.set_ylabel("Entries")
            ax.grid(True, alpha=0.25)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", ncol=3)
        fig.suptitle(f"{period}: elastic W-fit reach by electron angle")
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        fig.savefig(plot_dir / f"{safe}_theta_coverage.png", dpi=180)
        plt.close(fig)

        metrics = (
            ("peak_mean_error_gev", r"$\sigma(\mu_W)$ (GeV)",
             args.max_mean_error),
            ("peak_significance", r"Peak amplitude significance",
             args.minimum_significance),
            ("chi2_ndf", r"$\chi^2/\mathrm{ndf}$",
             args.max_chi2_ndf),
        )
        for column, ylabel, threshold in metrics:
            fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True)
            axes = axes.ravel()
            for sector in range(1, 7):
                ax = axes[sector - 1]
                sdf = period_df[period_df["sector"] == sector]
                ax.plot(sdf["theta_center_deg"], sdf[column], marker="o")
                ax.axhline(threshold, linestyle="--")
                ax.set_title(f"Sector {sector}")
                ax.set_xlabel(r"$\theta_e$ (deg)")
                ax.set_ylabel(ylabel)
                ax.grid(True, alpha=0.25)
            fig.suptitle(f"{period}: {ylabel} versus electron angle")
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            fig.savefig(plot_dir / f"{safe}_{column}.png", dpi=180)
            plt.close(fig)


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--su22-file", type=Path, default=DEFAULT_INPUTS["su22"])
    parser.add_argument("--fa22-file", type=Path, default=DEFAULT_INPUTS["fa22"])
    parser.add_argument("--sp23-file", type=Path, default=DEFAULT_INPUTS["sp23"])
    parser.add_argument("--tree-name", default="PhysicsEvents")
    parser.add_argument("--output-dir", type=Path,
                        default=Path("output/electron_theta_reach_scan"))
    parser.add_argument("--step-size", default="250 MB")
    parser.add_argument("--workers", type=int, default=4)

    parser.add_argument("--theta-min", type=float, default=7.0)
    parser.add_argument("--theta-max", type=float, default=25.0)
    parser.add_argument("--theta-bin-width", type=float, default=0.5)

    parser.add_argument("--w-preselection-min", type=float, default=0.65)
    parser.add_argument("--w-preselection-max", type=float, default=1.45)
    parser.add_argument("--w-hist-min", type=float, default=0.70)
    parser.add_argument("--w-hist-max", type=float, default=1.35)
    parser.add_argument("--w-fit-min", type=float, default=0.82)
    parser.add_argument("--w-fit-max", type=float, default=1.12)
    parser.add_argument("--w-bins", type=int, default=60)

    parser.add_argument("--minimum-events", type=int, default=500)
    parser.add_argument("--max-mean-error", type=float, default=0.010)
    parser.add_argument("--minimum-significance", type=float, default=3.0)
    parser.add_argument("--max-chi2-ndf", type=float, default=5.0)
    parser.add_argument("--min-sigma", type=float, default=0.008)
    parser.add_argument("--max-sigma", type=float, default=0.080)
    parser.add_argument("--min-peak-mean", type=float, default=0.86)
    parser.add_argument("--max-peak-mean", type=float, default=1.02)
    return parser


def main():
    args = build_parser().parse_args()

    inputs = {
        "su22": args.su22_file,
        "fa22": args.fa22_file,
        "sp23": args.sp23_file,
    }
    for path in inputs.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = args.output_dir / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    with ProcessPoolExecutor(max_workers=max(1, min(args.workers, 4))) as pool:
        future_map = {
            pool.submit(scan_period, period, inputs, args): period
            for period in PERIODS
        }
        for future in as_completed(future_map):
            period = future_map[future]
            print(f"[theta reach] Completed {period.label}")
            frames.append(future.result())

    results = classify(pd.concat(frames, ignore_index=True), args)
    results.to_csv(csv_dir / "theta_bin_fit_scan.csv", index=False)

    reach_records = []
    for (period, sector), group in results.groupby(["period", "sector"], sort=False):
        reach_records.append({
            "period": period,
            "sector": int(sector),
            "highest_contiguous_usable_theta_deg": contiguous_reach(group),
            "highest_any_usable_theta_deg": (
                float(group.loc[group["usable"], "theta_high_deg"].max())
                if group["usable"].any() else math.nan
            ),
            "usable_bins": int(group["usable"].sum()),
            "total_bins": int(len(group)),
        })

    reach = pd.DataFrame(reach_records)
    reach.to_csv(csv_dir / "theta_reach_summary.csv", index=False)

    make_plots(results, args.output_dir, args)

    print("\nHighest contiguous usable theta by period and sector")
    print(reach.to_string(index=False))
    print(f"\nOutputs written under {args.output_dir}")


if __name__ == "__main__":
    main()
