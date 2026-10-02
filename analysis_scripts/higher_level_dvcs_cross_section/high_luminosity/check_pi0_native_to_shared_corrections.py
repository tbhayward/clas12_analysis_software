#!/usr/bin/env python3

import numpy as np
import pandas as pd


INPUT = "output/pi0_gk_stage3/partons_gk/production_grid/gk_pi0_native_to_shared_corrections.csv"


def main():
    df = pd.read_csv(INPUT)

    print("Columns:")
    print(df.columns.tolist())
    print()

    required = {
        "reduced_value_native",
        "reduced_value_shared",
        "gk_shared_over_native",
    }
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(
            f"Correction table is missing reduced-cross-section columns: {sorted(missing)}"
        )
    #endif

    valid = (
        np.isfinite(df["reduced_value_native"])
        & np.isfinite(df["reduced_value_shared"])
        & np.isfinite(df["gk_shared_over_native"])
        & (df["reduced_value_native"] != 0.0)
    )
    direct = (
        df.loc[valid, "reduced_value_shared"].to_numpy(float)
        / df.loc[valid, "reduced_value_native"].to_numpy(float)
    )
    stored = df.loc[valid, "gk_shared_over_native"].to_numpy(float)
    max_residual = np.max(np.abs(stored / direct - 1.0)) if len(direct) else np.nan
    print(
        "Reduced-correction identity check: "
        f"max |stored/(reduced_shared/reduced_native)-1| = {max_residual:.3e}"
    )
    if not np.isfinite(max_residual) or max_residual > 1e-12:
        raise RuntimeError("Stored correction does not match the reduced PARTONS ratio.")
    #endif
    print()

    for campaign in ["rga", "rgk"]:
        x = df.loc[
            df["campaign"].str.lower() == campaign,
            "gk_shared_over_native",
        ].dropna().to_numpy(float)

        q = np.percentile(x, [0, 2.5, 16, 50, 84, 97.5, 100])

        print(campaign.upper())
        print(f"  N                  = {len(x)}")
        print(f"  median C           = {q[3]:.6f}")
        print(f"  central 68%        = [{q[2]:.6f}, {q[4]:.6f}]")
        print(f"  central 95%        = [{q[1]:.6f}, {q[5]:.6f}]")
        print(f"  min / max          = [{q[0]:.6f}, {q[6]:.6f}]")
        print(f"  median |C-1|       = {np.median(np.abs(x - 1)) * 100:.2f}%")
        print(f"  95th pct |C-1|     = {np.percentile(np.abs(x - 1), 95) * 100:.2f}%")
        print(f"  max |C-1|          = {np.max(np.abs(x - 1)) * 100:.2f}%")
        print(f"  fraction |C-1|> 5% = {np.mean(np.abs(x - 1) > 0.05) * 100:.1f}%")
        print(f"  fraction |C-1|>10% = {np.mean(np.abs(x - 1) > 0.10) * 100:.1f}%")
        print(f"  fraction |C-1|>20% = {np.mean(np.abs(x - 1) > 0.20) * 100:.1f}%")
        print()
    #endfor

    print("10 largest corrections away from unity:")
    z = df.dropna(subset=["gk_shared_over_native"]).copy()
    z["abs_percent_correction"] = 100 * np.abs(z["gk_shared_over_native"] - 1)
    cols = [
        "point_id",
        "campaign",
        "phi_deg",
        "gk_shared_over_native",
        "abs_percent_correction",
    ]
    print(z.nlargest(10, "abs_percent_correction")[cols].to_string(index=False))


if __name__ == "__main__":
    main()
#endif
