#!/usr/bin/env python3

import numpy as np
import pandas as pd
from pathlib import Path


CORRECTIONS = Path(
    "output/pi0_gk_stage3/partons_gk/production_grid/"
    "gk_pi0_native_to_shared_corrections.csv"
)

RGA = Path(
    "import/fa18_rosenbluth_inputs_20260924T165948Z/"
    "rga_10604/combined_reduced_cross_sections.csv"
)

RGK = Path(
    "import/fa18_rosenbluth_inputs_20260924T165948Z/"
    "rgk_6535/rgk6535_reduced_cross_sections.csv"
)

OUT = Path(
    "output/pi0_gk_stage3/diagnostics/"
    "actual_measured_phi_native_to_shared_corrections.csv"
)


def reconstruct_factors(cc, phi_deg):
    cc = cc.sort_values("phi_deg").copy()
    sample_phi = cc["phi_deg"].to_numpy(float)
    expected = np.array([0.0, 90.0, 180.0])

    if len(cc) != 3 or not np.allclose(sample_phi, expected, rtol=0.0, atol=1e-9):
        raise RuntimeError(
            f"Expected correction samples at 0, 90, 180 deg; got {sample_phi.tolist()}"
        )
    #endif

    p = np.deg2rad(sample_phi)
    H = np.column_stack([np.ones(3), np.cos(p), np.cos(2.0 * p)])

    native = np.linalg.solve(
        H, cc["partons_value_native"].to_numpy(float)
    )
    shared = np.linalg.solve(
        H, cc["partons_value_shared"].to_numpy(float)
    )

    pm = np.deg2rad(np.asarray(phi_deg, float))
    Hm = np.column_stack(
        [np.ones(len(pm)), np.cos(pm), np.cos(2.0 * pm)]
    )

    native_eval = Hm @ native
    shared_eval = Hm @ shared

    # Use a relative test based on the scale of this PARTONS curve itself.
    # Do not compare to an absolute O(1) scale: PARTONS observable values can
    # legitimately be much smaller than 1 in their reported units.
    scale = float(np.max(np.abs(native_eval)))
    if not np.isfinite(scale) or scale == 0.0:
        raise RuntimeError("Reconstructed native GK cross section is identically zero.")
    #endif

    tiny = np.finfo(float).eps * 100.0 * scale
    if np.any(np.abs(native_eval) <= tiny):
        raise RuntimeError(
            "Reconstructed native GK cross section is numerically zero "
            "relative to its own scale."
        )
    #endif

    factors = shared_eval / native_eval
    if np.any(~np.isfinite(factors)) or np.any(factors <= 0):
        raise RuntimeError("Invalid reconstructed native-to-shared factor.")
    #endif

    return factors


def summarize(x, label):
    x = np.asarray(x, float)
    a = np.abs(x - 1.0)
    q = np.percentile(x, [0, 2.5, 16, 50, 84, 97.5, 100])

    print(label)
    print(f"  N                   = {len(x)}")
    print(f"  median C            = {q[3]:.6f}")
    print(f"  central 68%         = [{q[2]:.6f}, {q[4]:.6f}]")
    print(f"  central 95%         = [{q[1]:.6f}, {q[5]:.6f}]")
    print(f"  min / max           = [{q[0]:.6f}, {q[6]:.6f}]")
    print(f"  median |C-1|        = {np.median(a) * 100:.2f}%")
    print(f"  95th pct |C-1|      = {np.percentile(a, 95) * 100:.2f}%")
    print(f"  max |C-1|           = {np.max(a) * 100:.2f}%")
    print(f"  fraction |C-1| > 5% = {np.mean(a > 0.05) * 100:.1f}%")
    print(f"  fraction |C-1|>10%  = {np.mean(a > 0.10) * 100:.1f}%")
    print(f"  fraction |C-1|>20%  = {np.mean(a > 0.20) * 100:.1f}%")
    print()


def main():
    corr = pd.read_csv(CORRECTIONS)
    rows = []

    for campaign, path in [("rga", RGA), ("rgk", RGK)]:
        data = pd.read_csv(path)

        required = {"iq2", "ixb", "it", "iphi", "phi_center_deg"}
        missing = required - set(data.columns)
        if missing:
            raise RuntimeError(
                f"{campaign.upper()} input is missing columns: {sorted(missing)}"
            )
        #endif

        # The R#### point IDs are defined by the same sorted common-cell ordering
        # used by the Stage-1/Stage-3 pipeline.
        cells = (
            data[["iq2", "ixb", "it"]]
            .drop_duplicates()
            .sort_values(["iq2", "ixb", "it"])
        )

        # Build the common-cell ordering from both campaigns below, rather than
        # assuming either campaign alone contains every common cell.
        if campaign == "rga":
            rga_data = data
        else:
            rgk_data = data
        #endif
    #endfor

    rga_cells = rga_data[["iq2", "ixb", "it"]].drop_duplicates()
    rgk_cells = rgk_data[["iq2", "ixb", "it"]].drop_duplicates()
    common = (
        rga_cells.merge(rgk_cells, on=["iq2", "ixb", "it"], how="inner")
        .sort_values(["iq2", "ixb", "it"])
        .reset_index(drop=True)
    )
    common["point_id"] = [f"R{i:04d}" for i in range(len(common))]

    for campaign, data in [("rga", rga_data), ("rgk", rgk_data)]:
        for r in common.itertuples(index=False):
            pid = r.point_id
            cc = corr[
                (corr["point_id"].astype(str) == pid)
                & (corr["campaign"].astype(str).str.lower() == campaign)
            ].copy()

            # Points excluded before/while producing the correction grid simply
            # have no correction rows and are not part of this diagnostic.
            if cc.empty:
                continue
            #endif

            g = data[
                (data["iq2"] == r.iq2)
                & (data["ixb"] == r.ixb)
                & (data["it"] == r.it)
            ].copy()

            if g.empty:
                continue
            #endif

            try:
                fac = reconstruct_factors(
                    cc, g["phi_center_deg"].to_numpy(float)
                )
            except RuntimeError as exc:
                print(
                    f"SKIP {pid} {campaign.upper()} "
                    f"(iq2={r.iq2}, ixb={r.ixb}, it={r.it}): {exc}"
                )
                continue
            #endtry

            for (_, drow), f in zip(g.iterrows(), fac):
                rows.append(
                    {
                        "point_id": pid,
                        "campaign": campaign,
                        "iq2": int(r.iq2),
                        "ixb": int(r.ixb),
                        "it": int(r.it),
                        "iphi": int(drow["iphi"]),
                        "phi_deg": float(drow["phi_center_deg"]),
                        "gk_shared_over_native": float(f),
                        "fractional_shift": float(f - 1.0),
                        "abs_percent_shift": float(100.0 * abs(f - 1.0)),
                    }
                )
            #endfor
        #endfor
    #endfor

    out = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)

    print(f"Wrote {OUT}")
    print(f"Total measured phi points with corrections: {len(out)}")
    print()

    for campaign in ["rga", "rgk"]:
        x = out.loc[
            out["campaign"] == campaign, "gk_shared_over_native"
        ].to_numpy(float)
        summarize(x, campaign.upper())
    #endfor

    print("Correction size by measured phi region:")
    bins = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360]
    out["phi_region"] = pd.cut(
        np.mod(out["phi_deg"], 360.0),
        bins=bins,
        include_lowest=True,
        right=False,
    )

    table = (
        out.groupby(["campaign", "phi_region"], observed=True)
        .agg(
            N=("gk_shared_over_native", "size"),
            median_C=("gk_shared_over_native", "median"),
            median_abs_percent=("abs_percent_shift", "median"),
            max_abs_percent=("abs_percent_shift", "max"),
        )
    )
    print(table.to_string())
    print()

    print("10 largest corrections actually applied to measured phi bins:")
    cols = [
        "point_id",
        "campaign",
        "iphi",
        "phi_deg",
        "gk_shared_over_native",
        "abs_percent_shift",
    ]
    print(
        out.nlargest(10, "abs_percent_shift")[cols].to_string(index=False)
    )
    print()

    print("10 cells with largest maximum correction:")
    cell = (
        out.groupby(["point_id", "campaign"], as_index=False)
        .agg(
            n_phi=("gk_shared_over_native", "size"),
            median_C=("gk_shared_over_native", "median"),
            median_abs_percent=("abs_percent_shift", "median"),
            max_abs_percent=("abs_percent_shift", "max"),
        )
        .sort_values("max_abs_percent", ascending=False)
    )
    print(cell.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
