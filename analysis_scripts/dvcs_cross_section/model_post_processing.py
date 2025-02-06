#!/usr/bin/env python3

import sys
import os
import subprocess
import numpy as np
import pandas as pd

# GEPARD imports
import gepard as g
from gepard.fits import th_KM15

###############################################
# 0) Command-Line Argument for Number of Rows
###############################################

def parse_args():
    """
    Parses sys.argv for an optional integer specifying how many rows to process.
    If none provided, returns None (meaning process all).
    """
    if len(sys.argv) > 1:
        try:
            return int(sys.argv[1])
        except ValueError:
            pass
    return None

###############################################
# 1) Define helper functions for KM15 & dvcsgen
###############################################

def km15_model(xB, Q2, t_pos, phi_deg, beam_E=10.604):
    """
    Compute KM15 cross section (via GEPARD) for given kinematics.

    - xB, Q2   : usual DIS variables
    - t_pos    : positive number for |t| (we feed negative to KM15 internally)
    - phi_deg  : phi in degrees, converted to radians with (pi - phi_rad) for TRENTO
    - beam_E   : beam energy (GeV), e.g. 10.604 or 10.1998
    """
    t_km15 = -abs(t_pos)  # Convert to negative t for KM15
    phi_rad = np.radians(phi_deg)
    phi_trento = np.pi - phi_rad  # TRENTO shift

    pt = g.DataPoint(
        xB        = xB,
        t         = t_km15,
        Q2        = Q2,
        phi       = phi_trento,
        observable= 'XS',
        frame     = 'trento',
        process   = 'ep2epgamma',
        exptype   = 'fixed target',
        in1energy = beam_E,
        in1charge = -1,
        in1polarization = 0
    )
    pt.prepare()
    return th_KM15.predict(pt)
#endfor


def dvcsgen_vgg(xB, Q2, t_pos, phi_deg, beam_E=10.604, globalfit=True, pol=0, local=False):
    """
    Calls external dvcsgen for the VGG model.
    """
    my_env = os.environ.copy()
    path = "/u/home/thayward/dvcsgen"
    if local:
        path = "/Users/sangbaek.lee/CLAS12/dvcs/print/"

    my_env["PATH"] = f"{path}:{my_env['PATH']}"
    my_env["CLASDVCS_PDF"] = path

    phi_rad = np.radians(phi_deg)
    cmd = [
        f"{path}/dvcsgen",
        "--beam", f"{beam_E:.3f}",
        "--x", f"{xB}", f"{xB}",
        "--q2", f"{Q2}", f"{Q2}",
        "--t", f"{t_pos}", f"{t_pos}",
        "--bh", "3",
        "--phi", f"{phi_rad}",
        "--gpd", "101",
        "--ycol", "0.0001"
    ]
    if globalfit:
        cmd.append("--globalfit")

    try:
        dstot = subprocess.check_output(cmd, env=my_env)
        # pol=0 => last line index=0; pol=1 => index=2; etc.
        if pol == 0:
            i = 0
        elif pol == 1:
            i = 2
        else:
            i = 1

        val_str = dstot.splitlines()[-1 - i].decode("utf-8")
        return float(val_str)
    except Exception as e:
        print(f"dvcsgen VGG error at xB={xB}, Q2={Q2}, t={t_pos}, phi={phi_deg} deg -> {e}")
        return 0.0
#endfor


def dvcsgen_bh_only(xB, Q2, t_pos, phi_deg, beam_E=10.604, globalfit=True, local=False):
    """
    Calls external dvcsgen for BH-only (bh=1).
    """
    my_env = os.environ.copy()
    path = "/u/home/thayward/dvcsgen"
    if local:
        path = "/Users/sangbaek.lee/CLAS12/dvcs/print/"

    my_env["PATH"] = f"{path}:{my_env['PATH']}"
    my_env["CLASDVCS_PDF"] = path

    phi_rad = np.radians(phi_deg)
    cmd = [
        f"{path}/dvcsgen",
        "--beam", f"{beam_E:.3f}",
        "--x", f"{xB}", f"{xB}",
        "--q2", f"{Q2}", f"{Q2}",
        "--t", f"{t_pos}", f"{t_pos}",
        "--bh", "1",
        "--phi", f"{phi_rad}",
        "--ycol", "0.0001"
    ]
    if globalfit:
        cmd.append("--globalfit")

    try:
        dstot = subprocess.check_output(cmd, env=my_env)
        val_str = dstot.splitlines()[-1].decode("utf-8")
        return float(val_str)
    except Exception as e:
        print(f"dvcsgen BH-only error at xB={xB}, Q2={Q2}, t={t_pos}, phi={phi_deg} deg -> {e}")
        return 0.0
#endfor

###############################################
# 1B) dvcsgen with --printrad to get Frad
###############################################

def dvcsgen_printrad(xB, Q2, t_pos, phi_deg, beam_E=10.604):
    """
    Calls dvcsgen with a special command:
       --printrad
    and parses the penultimate line for 'Frad_with_error' factor + sys.
    We then return (Frad_factor, Frad_sys).
    If something fails or doesn't parse, return (1.0, 0.0).

    Now includes debug print statements to show all dvcsgen output.
    """
    my_env = os.environ.copy()
    path = "/u/home/thayward/dvcsgen"

    my_env["PATH"] = f"{path}:{my_env['PATH']}"
    my_env["CLASDVCS_PDF"] = path

    phi_rad = np.radians(phi_deg)

    cmd = [
        f"{path}/dvcsgen",
        "--beam", f"{beam_E:.3f}",
        "--x", str(xB), str(xB),
        "--q2", str(Q2), str(Q2),
        "--t", str(t_pos), str(t_pos),
        "--gpd", "101",
        "--y", "0", "1",
        "--phi", f"{phi_rad}",
        "--vv2cut", "0.3",
        "--delta", "0.1",
        "--printrad"
    ]

    try:
        dstot = subprocess.check_output(cmd, env=my_env)
        decoded = dstot.decode("utf-8", errors="replace")

        # DEBUG PRINT: Full dvcsgen output
        print("===== dvcsgen printrad DEBUG OUTPUT =====")
        print(f"xB={xB}, Q2={Q2}, t={t_pos}, phi_deg={phi_deg}, beamE={beam_E}")
        print("CMD =", " ".join(cmd))
        print("---- dvcsgen output ----")
        print(decoded)
        print("===== END dvcsgen printrad DEBUG =====")

        lines = decoded.splitlines()
        if len(lines) < 2:
            print(" -> Too few lines in dvcsgen output, returning fallback (1.0, 0.0).")
            return (1.0, 0.0)

        # The penultimate line is lines[-2].
        # It should contain something like:
        # " Frad_with_error  0.88183283631317655  7.6579640022232884E-005"
        penultimate = lines[-2].strip()
        if "Frad_with_error" not in penultimate:
            print(f" -> 'Frad_with_error' not found in penultimate line: {penultimate}")
            return (1.0, 0.0)

        tokens = penultimate.split()
        # e.g. ["Frad_with_error","0.88183283631317655","7.6579640022232884E-005"]
        if len(tokens) < 3:
            print(f" -> penultimate line does not have 3 tokens: {penultimate}")
            return (1.0, 0.0)

        factor = float(tokens[1])
        sysval = float(tokens[2])
        print(f" -> PARSED Frad factor={factor:.6f}, sys={sysval:.6f}")
        return (factor, sysval)

    except subprocess.CalledProcessError as e:
        # If dvcsgen fails to run, we can see the error code or partial output
        print(f"dvcsgen printrad error (xB={xB}, Q2={Q2}, t={t_pos}, phi={phi_deg}): {e}")
        return (1.0, 0.0)
    except Exception as e:
        # Catch-all for any other parse/time out/etc.
        print(f"dvcsgen printrad unknown error (xB={xB}, Q2={Q2}, t={t_pos}, phi={phi_deg}): {e}")
        return (1.0, 0.0)
#enddef

###############################################
# 2) Function for bin-centering correction (KM15 + VGG)
###############################################
def calculate_fbin(row, prefix, beam_E, n_steps=3):
    """
    Same as your existing Fbin function: sub-binning with KM15 & VGG.
    Returns: (km15_fbin, vgg_fbin, final_fbin, fbin_sys_unc).
    """
    Mp = 0.938272
    xB_samples   = np.linspace(row['xB_min'],   row['xB_max'],   n_steps)
    Q2_samples   = np.linspace(row['Q2_min'],   row['Q2_max'],   n_steps)
    t_pos_samples= np.linspace(row['t_min'],    row['t_max'],    n_steps)
    phi_samples  = np.linspace(row['phi_min'],  row['phi_max'],  n_steps)

    valid_KM15 = []
    valid_VGG  = []

    for xB in xB_samples:
        for Q2 in Q2_samples:
            for t_pos in t_pos_samples:
                t_phys = -abs(t_pos)
                # t_min
                try:
                    sqrt_term = np.sqrt(1 + (4 * Mp**2 * xB**2) / Q2)
                    t_min_val = -Q2 * (1 - xB)**2 / (xB * (1 + sqrt_term))
                except:
                    continue

                try:
                    y = Q2 / (2 * Mp * xB * beam_E)
                    W = np.sqrt(Mp**2 + Q2 * (1/xB - 1))
                except:
                    continue

                if (t_phys >= t_min_val) and (0.19 < y < 0.8) and (W > 2.0):
                    for phi_deg in phi_samples:
                        try:
                            km15_val = km15_model(xB, Q2, t_pos, phi_deg, beam_E)
                            vgg_val  = dvcsgen_vgg(xB, Q2, t_pos, phi_deg, beam_E, globalfit=False)
                            valid_KM15.append(km15_val)
                            valid_VGG.append(vgg_val)
                        except:
                            continue
                else:
                    continue

    # Center values
    try:
        centerKM15 = km15_model(row['xB_avg'], row['Q2_avg'], row['t_avg'], row['phi_avg'], beam_E)
        centerVGG  = dvcsgen_vgg(row['xB_avg'], row['Q2_avg'], row['t_avg'], row['phi_avg'], beam_E, globalfit=False)
    except:
        return (1.0, 1.0, 1.0, 0.0)

    if not valid_KM15 or not valid_VGG:
        return (1.0, 1.0, 1.0, 0.0)

    avgKM15 = np.mean(valid_KM15)
    avgVGG  = np.mean(valid_VGG)
    if avgKM15 == 0.0 or avgVGG == 0.0:
        return (1.0, 1.0, 1.0, 0.0)

    km15_fbin = centerKM15 / avgKM15
    vgg_fbin  = centerVGG  / avgVGG
    fbin_values   = [km15_fbin, vgg_fbin]
    final_fbin    = np.mean(fbin_values)
    fbin_sys_unc  = np.std(fbin_values)

    return (km15_fbin, vgg_fbin, final_fbin, fbin_sys_unc)
#enddef

###############################################
# 2B) Radiative correction (Frad) with sub-binning
###############################################
def calculate_frad(row, prefix, beam_E, n_steps=3):
    """
    Similar sub-binning approach, but calls dvcsgen with `--printrad`.

    Procedure:
      - For each sub-bin, call dvcsgen_printrad(...), gather sub-bin factor (f_sub).
      - Then at bin center, also call dvcsgen_printrad(...) => (f_center, sys_center).
      - final_frad = average_subbin / f_center
      - final_frad_sys = sqrt( (std_subbin^2) + (sys_center^2 ) )

    returns (final_frad, final_frad_sys).
    """
    Mp = 0.938272

    xB_samples   = np.linspace(row['xB_min'],   row['xB_max'],   n_steps)
    Q2_samples   = np.linspace(row['Q2_min'],   row['Q2_max'],   n_steps)
    t_pos_samples= np.linspace(row['t_min'],    row['t_max'],    n_steps)
    phi_samples  = np.linspace(row['phi_min'],  row['phi_max'],  n_steps)

    subbin_vals = []  # store Frad for each sub-bin

    for xB in xB_samples:
        for Q2 in Q2_samples:
            for t_pos in t_pos_samples:
                t_phys = -abs(t_pos)
                try:
                    sqrt_term = np.sqrt(1 + (4 * Mp**2 * xB**2) / Q2)
                    t_min_val = -Q2 * (1 - xB)**2 / (xB * (1 + sqrt_term))
                except:
                    continue

                try:
                    y = Q2 / (2*Mp*xB*beam_E)
                    W = np.sqrt(Mp**2 + Q2*(1/xB - 1))
                except:
                    continue

                if (t_phys >= t_min_val) and (0.19 < y < 0.8) and (W>2.0):
                    for phi_deg in phi_samples:
                        try:
                            f_sub, _subsys = dvcsgen_printrad(xB, Q2, t_pos, phi_deg, beam_E)
                            subbin_vals.append(f_sub)
                        except:
                            continue
                else:
                    continue

    # Now the center
    try:
        f_center, center_sys = dvcsgen_printrad(row['xB_avg'], row['Q2_avg'], row['t_avg'], row['phi_avg'], beam_E)
    except:
        return (1.0, 0.0)

    if (not subbin_vals) or (f_center==0.0):
        return (1.0, 0.0)

    avg_subbin = np.mean(subbin_vals)
    std_subbin = np.std(subbin_vals)

    # final factor = avg_subbin / f_center
    # total sys = sqrt( std_subbin^2 + center_sys^2 )
    final_val = avg_subbin / f_center
    final_sys = np.sqrt(std_subbin**2 + center_sys**2)

    return (final_val, final_sys)
#enddef

###################################################
# 3) Main code
###################################################
def main():
    # 3A) parse arguments
    n_rows_to_process = parse_args()  # None if no argument provided

    print("Beginning code.")
    input_csv = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/output/unfolding_data.csv"
    print(f"Reading CSV from: {input_csv}")
    df = pd.read_csv(input_csv)
    print(f"Loaded DataFrame with {len(df)} rows and {len(df.columns)} columns.")

    # If user specified a limit, override
    if n_rows_to_process is None:
        n_rows_to_process = len(df)
    else:
        n_rows_to_process = min(n_rows_to_process, len(df))

    # Create columns for bin-centering factors
    fb_cols = ['KM15_Fbin', 'VGG_Fbin', 'Fbin', 'Fbin_sys_uncertainty']
    # Create columns for radiative corrections
    frad_cols = ['Frad', 'Frad_sys_uncertainty']

    print("Preparing new columns for Fbin and Frad...")

    for prefix in ['fall', 'spring']:
        for c in fb_cols:
            df[f"{prefix}_{c}"] = np.nan
        for c in frad_cols:
            df[f"{prefix}_{c}"] = np.nan

    # Also ensure we have columns for existing model calculations
    # (km15, vgg, bh) if not present
    # Not strictly necessary, but in case they don't exist
    for prefix in ['fall', 'spring']:
        for c in ['km15','vgg','bh']:
            colname = f"{prefix}_{c}"
            if colname not in df.columns:
                df[colname] = np.nan

    print("Beginning cross-section, Fbin, and Frad calculations...")

    # You may choose a progress interval
    progress_interval = 25

    for i in range(n_rows_to_process):
        if i % progress_interval == 0:
            print(f"  Processing row {i} of {n_rows_to_process} ...")

        row = df.iloc[i]

        # 1) "Existing" model calculations at bin center
        fall_km15_val = km15_model(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.604)
        fall_vgg_val  = dvcsgen_vgg(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.604, globalfit=False)
        fall_bh_val   = dvcsgen_bh_only(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.604, globalfit=False)

        spring_km15_val = km15_model(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.1998)
        spring_vgg_val  = dvcsgen_vgg(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.1998, globalfit=False)
        spring_bh_val   = dvcsgen_bh_only(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.1998, globalfit=False)

        # 2) Bin-centering factor (Fall)
        fall_KM15_Fbin, fall_VGG_Fbin, fall_Fbin, fall_Fbin_sys = calculate_fbin(row, 'fall', 10.604)
        df.loc[i, "fall_KM15_Fbin"]            = fall_KM15_Fbin
        df.loc[i, "fall_VGG_Fbin"]             = fall_VGG_Fbin
        df.loc[i, "fall_Fbin"]                 = fall_Fbin
        df.loc[i, "fall_Fbin_sys_uncertainty"] = fall_Fbin_sys

        # 2B) Radiative correction factor (Fall)
        fall_Frad_val, fall_Frad_sys = calculate_frad(row, 'fall', 10.604)
        df.loc[i, "fall_Frad"] = fall_Frad_val
        df.loc[i, "fall_Frad_sys_uncertainty"] = fall_Frad_sys

        # 3) Bin-centering factor (Spring)
        spring_KM15_Fbin, spring_VGG_Fbin, spring_Fbin, spring_Fbin_sys = calculate_fbin(row, 'spring', 10.1998)
        df.loc[i, "spring_KM15_Fbin"]            = spring_KM15_Fbin
        df.loc[i, "spring_VGG_Fbin"]             = spring_VGG_Fbin
        df.loc[i, "spring_Fbin"]                 = spring_Fbin
        df.loc[i, "spring_Fbin_sys_uncertainty"] = spring_Fbin_sys

        # 3B) Radiative correction factor (Spring)
        spring_Frad_val, spring_Frad_sys = calculate_frad(row, 'spring', 10.1998)
        df.loc[i, "spring_Frad"] = spring_Frad_val
        df.loc[i, "spring_Frad_sys_uncertainty"] = spring_Frad_sys

        # 4) Store bin-center model outputs
        df.loc[i, "fall_km15"]   = fall_km15_val
        df.loc[i, "fall_vgg"]    = fall_vgg_val
        df.loc[i, "fall_bh"]     = fall_bh_val
        df.loc[i, "spring_km15"] = spring_km15_val
        df.loc[i, "spring_vgg"]  = spring_vgg_val
        df.loc[i, "spring_bh"]   = spring_bh_val

        # 5) Apply final Fbin *and* Frad factors to cross section & stat uncertainty.
        for pre in ['fall', 'spring']:
            fbin_factor = df.loc[i, f"{pre}_Fbin"]
            frad_factor = df.loc[i, f"{pre}_Frad"]
            
            total_factor = fbin_factor * frad_factor
            
            df.loc[i, f"{pre}_cross_section"] *= total_factor
            df.loc[i, f"{pre}_cross_section_stat_uncertainty"] *= total_factor

    # Reordering columns
    print("Reordering columns...")
    all_cols = list(df.columns)

    def get_ordered_columns(prefix):
        """
        Return columns in correct order for a given prefix.
        We'll place the new Frad columns right after Fbin_sys_uncertainty.
        """
        return [
            f"{prefix}_bin_volume",
            f"{prefix}_KM15_Fbin",
            f"{prefix}_VGG_Fbin",
            f"{prefix}_Fbin",
            f"{prefix}_Fbin_sys_uncertainty",
            f"{prefix}_Frad",
            f"{prefix}_Frad_sys_uncertainty",
            f"{prefix}_cross_section",
            f"{prefix}_cross_section_stat_uncertainty",
            f"{prefix}_cross_section_sys_uncertainty",
            f"{prefix}_km15",
            f"{prefix}_vgg",
            f"{prefix}_bh"
        ]

    new_column_order = []
    common_columns   = [c for c in all_cols if not c.startswith(("fall_","spring_"))]

    fall_columns   = get_ordered_columns('fall')
    spring_columns = get_ordered_columns('spring')

    # keep original non-prefixed first
    for col in common_columns:
        new_column_order.append(col)
        if col in all_cols:
            all_cols.remove(col)

    # then fall
    for col in fall_columns:
        if col in all_cols:
            new_column_order.append(col)
            all_cols.remove(col)

    # then spring
    for col in spring_columns:
        if col in all_cols:
            new_column_order.append(col)
            all_cols.remove(col)

    # any leftover
    new_column_order += all_cols

    df = df[new_column_order]

    # Save
    outpath = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/output/unfolding_data_with_models.csv"
    df.to_csv(outpath, index=False)
    print(f"Done! Updated file saved to:\n{outpath}")


###############################################
# 4) Entry point
###############################################
if __name__ == "__main__":
    main()