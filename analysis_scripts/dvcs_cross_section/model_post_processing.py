import pandas as pd
import numpy as np
import os
import subprocess

# GEPARD imports
import gepard as g
from gepard.fits import th_KM15

###############################################
# 1) Define helper functions for KM15 & dvcsgen
###############################################

def km15_model(xB, Q2, t_pos, phi_deg, beam_E=10.604):
    """
    Compute KM15 cross section (via GEPARD) for given kinematics.

    - xB, Q2   : usual DIS variables
    - t_pos    : positive number for |t|, we feed it as negative to KM15
    - phi_deg  : phi in degrees, converted to radians with (pi - phi_rad) for TRENTO
    - beam_E   : beam energy (GeV), e.g. 10.604 or 10.1998
    """
    # Convert to negative t for KM15
    t_km15 = -abs(t_pos)

    # Convert degrees to radians and shift for TRENTO frame
    phi_rad = np.radians(phi_deg)
    phi_trento = np.pi - phi_rad

    # Construct GEPARD DataPoint
    pt = g.DataPoint(
        xB=xB,
        t=t_km15,
        Q2=Q2,
        phi=phi_trento,
        observable='XS',
        frame='trento',
        process='ep2epgamma',
        exptype='fixed target',
        in1energy=beam_E,
        in1charge=-1,
        in1polarization=0
    )
    pt.prepare()
    return th_KM15.predict(pt)
#endfor

def dvcsgen_vgg(xB, Q2, t_pos, phi_deg, beam_E=10.604, globalfit=True, pol=0, local=False):
    """
    Calls an external dvcsgen executable for the VGG model (unused for bin-centering if desired,
    but here it *is* used to define sub-binning for a second model, as we want the final Fbin
    to be the average of KM15 & VGG).
    """
    my_env = os.environ.copy()
    path = "/u/home/thayward/dvcsgen"
    if local:
        path = "/Users/sangbaek.lee/CLAS12/dvcs/print/"

    my_env["PATH"] = f"{path}:{my_env['PATH']}"
    my_env["CLASDVCS_PDF"] = path

    # dvcsgen expects phi in radians
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
    Calls an external dvcsgen executable for the BH-only model (bh=1).
    (Unused for bin-centering, but kept if you want BH in existing calculations.)
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

###################################################
# 2) Function for bin-centering correction (KM15 + VGG)
###################################################
def calculate_fbin(row, prefix, beam_E, n_steps=5):
    """
    Calculate bin-centering factors by sub-binning with both KM15 and VGG.
    We return four values:

        (km15_fbin, vgg_fbin, final_fbin, fbin_sys_unc)

    Where:
      - km15_fbin = centerKM15 / average over sub-bins (KM15)
      - vgg_fbin  = centerVGG  / average over sub-bins (VGG)
      - final_fbin = (km15_fbin + vgg_fbin)/2
      - fbin_sys_unc = std([km15_fbin, vgg_fbin])

    This final_fbin is what we'll store in row["prefix_Fbin"] and apply to the cross section.
    """
    Mp = 0.938272  # Proton mass in GeV/c²

    # Generate grid points for each variable
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

                # y, W
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

    # Now compute the center values (KM15 & VGG) at bin center
    try:
        centerKM15 = km15_model(row['xB_avg'], row['Q2_avg'], row['t_avg'], row['phi_avg'], beam_E)
        centerVGG  = dvcsgen_vgg(row['xB_avg'], row['Q2_avg'], row['t_avg'], row['phi_avg'], beam_E, globalfit=False)
    except:
        # fallback
        return (1.0, 1.0, 1.0, 0.0)

    if not valid_KM15 or not valid_VGG:
        # no valid sub-bins
        return (1.0, 1.0, 1.0, 0.0)

    avgKM15 = np.mean(valid_KM15)
    avgVGG  = np.mean(valid_VGG)
    if avgKM15 == 0.0 or avgVGG == 0.0:
        return (1.0, 1.0, 1.0, 0.0)

    # compute each ratio
    km15_fbin = centerKM15 / avgKM15
    vgg_fbin  = centerVGG  / avgVGG
    # final fbin is average, sys is the std
    fbin_values   = [km15_fbin, vgg_fbin]
    final_fbin    = np.mean(fbin_values)
    fbin_sys_unc  = np.std(fbin_values)

    # debug prints if desired:
    print(f"KM15: {km15_fbin:.4f}, VGG: {vgg_fbin:.4f}, final: {final_fbin:.4f}, sys: {fbin_sys_unc:.4f}")

    return (km15_fbin, vgg_fbin, final_fbin, fbin_sys_unc)
#enddef

###################################################
# 3) Main code execution
###################################################
if __name__ == "__main__":
    print("Beginning code.")
    input_csv = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/output/unfolding_data.csv"
    print(f"Reading CSV from: {input_csv}")
    df = pd.read_csv(input_csv)
    print(f"Loaded DataFrame with {len(df)} rows and {len(df.columns)} columns.")

    # Create columns for the Fbin factors (KM15_Fbin, VGG_Fbin, Fbin, Fbin_sys_uncertainty)
    print("Preparing new columns...")

    fb_cols = ['KM15_Fbin', 'VGG_Fbin', 'Fbin', 'Fbin_sys_uncertainty']
    for prefix in ['fall', 'spring']:
        for c in fb_cols:
            df[f"{prefix}_{c}"] = np.nan

    # Main loop
    print("Beginning cross-section and Fbin calculations...")
    progress_interval = 1  # or whatever you prefer

    for i in range(len(df)):
        if i % progress_interval == 0:
            print(f"  Processing row {i} of {len(df)} ...")

        row = df.iloc[i]

        # 1) "Existing" model calculations, if you want them:
        #    (These are the cross section predictions at the bin center.)
        fall_km15_val = km15_model(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.604)
        fall_vgg_val  = dvcsgen_vgg(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.604, globalfit=False)
        fall_bh_val   = dvcsgen_bh_only(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.604, globalfit=False)

        spring_km15_val = km15_model(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.1998)
        spring_vgg_val  = dvcsgen_vgg(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.1998, globalfit=False)
        spring_bh_val   = dvcsgen_bh_only(row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.1998, globalfit=False)

        # 2) Bin-centering factor (Fall)
        fall_KM15_Fbin, fall_VGG_Fbin, fall_Fbin, fall_Fbin_sys = calculate_fbin(row, 'fall', beam_E=10.604)
        # store them
        df.loc[i, "fall_KM15_Fbin"]              = fall_KM15_Fbin
        df.loc[i, "fall_VGG_Fbin"]               = fall_VGG_Fbin
        df.loc[i, "fall_Fbin"]                   = fall_Fbin
        df.loc[i, "fall_Fbin_sys_uncertainty"]   = fall_Fbin_sys

        # 3) Bin-centering factor (Spring)
        spring_KM15_Fbin, spring_VGG_Fbin, spring_Fbin, spring_Fbin_sys = calculate_fbin(row, 'spring', beam_E=10.1998)
        df.loc[i, "spring_KM15_Fbin"]            = spring_KM15_Fbin
        df.loc[i, "spring_VGG_Fbin"]             = spring_VGG_Fbin
        df.loc[i, "spring_Fbin"]                 = spring_Fbin
        df.loc[i, "spring_Fbin_sys_uncertainty"] = spring_Fbin_sys

        # 4) Store the bin-center model outputs
        df.loc[i, "fall_km15"]   = fall_km15_val
        df.loc[i, "fall_vgg"]    = fall_vgg_val
        df.loc[i, "fall_bh"]     = fall_bh_val
        df.loc[i, "spring_km15"] = spring_km15_val
        df.loc[i, "spring_vgg"]  = spring_vgg_val
        df.loc[i, "spring_bh"]   = spring_bh_val

        # 5) Apply final Fbin to cross section & stat. uncertainty
        #    (We use 'fall_Fbin' or 'spring_Fbin' = the average of KM15 & VGG.)
        for pre in ['fall', 'spring']:
            fbin_val = df.loc[i, f"{pre}_Fbin"]
            df.loc[i, f"{pre}_cross_section"] *= fbin_val
            df.loc[i, f"{pre}_cross_section_stat_uncertainty"] *= fbin_val

    ###################################################
    # Reorder columns
    ###################################################
    print("Reordering columns...")
    all_cols = list(df.columns)

    def get_ordered_columns(prefix):
        """Return columns in correct order for a given prefix."""
        return [
            f"{prefix}_bin_volume",
            f"{prefix}_KM15_Fbin",
            f"{prefix}_VGG_Fbin",
            f"{prefix}_Fbin",
            f"{prefix}_Fbin_sys_uncertainty",
            f"{prefix}_cross_section",
            f"{prefix}_cross_section_stat_uncertainty",
            f"{prefix}_cross_section_sys_uncertainty",
            f"{prefix}_km15",
            f"{prefix}_vgg",
            f"{prefix}_bh"
        ]

    new_column_order = []
    common_columns   = [col for col in all_cols if not col.startswith(('fall_', 'spring_'))]

    # Add columns for FALL, then SPRING, in desired order
    fall_columns   = get_ordered_columns('fall')
    spring_columns = get_ordered_columns('spring')

    # Preserve original non-prefixed columns at start
    for col in common_columns:
        new_column_order.append(col)
        if col in all_cols:
            all_cols.remove(col)

    # Then add fall columns
    for col in fall_columns:
        if col in all_cols:
            new_column_order.append(col)
            all_cols.remove(col)

    # Then add spring columns
    for col in spring_columns:
        if col in all_cols:
            new_column_order.append(col)
            all_cols.remove(col)

    # Add any leftover
    new_column_order += all_cols

    # Reorder
    df = df[new_column_order]

    # Save output
    output_csv = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/output/unfolding_data_with_models.csv"
    df.to_csv(output_csv, index=False)
    print(f"Done! Updated file saved to:\n{output_csv}")