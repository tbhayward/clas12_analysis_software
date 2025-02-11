#!/usr/bin/env python3

import sys
import os
import subprocess
import numpy as np
import pandas as pd
import concurrent.futures  # NEW for parallel sub-binning

# GEPARD imports
import gepard as g
from gepard.fits import th_KM15

######################################################
# 0) Command-Line Argument for Number of Rows
######################################################

def parse_args():
    """
    Parses sys.argv for an optional integer specifying how many rows to process.
    If none provided, returns None (meaning process all).
    Example:
      python model_post_processing.py 100
    will process only the first 100 rows.
    """
    if len(sys.argv) > 1:
        try:
            return int(sys.argv[1])
        except ValueError:
            pass
    return None

######################################################
# PATHS FOR SEPARATE DVCSGEN INSTALLATIONS
######################################################
REGULAR_DVCSGEN_PATH = "/u/home/thayward/dvcsgens/dvcsgen_print"
RADIATIVE_DVCSGEN_PATH = "/u/home/thayward/dvcsgens/dvcsgen_rad2"

######################################################
# 1) Define helper functions for KM15 & dvcsgen
######################################################

def km15_model(xB, Q2, t_pos, phi_deg, beam_E=10.604):
    t_km15 = -abs(t_pos)
    phi_rad = np.radians(phi_deg)
    phi_trento = np.pi - phi_rad

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
    my_env = os.environ.copy()
    path = REGULAR_DVCSGEN_PATH
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
        if pol == 0:
            i = 0
        elif pol == 1:
            i = 2
        else:
            i = 1

        line_str = dstot.splitlines()[-1 - i].decode("utf-8")
        tokens = line_str.split()
        val_str = tokens[-1]
        return float(val_str)
    except Exception as e:
        print(f"dvcsgen VGG error at xB={xB}, Q2={Q2}, t={t_pos}, phi={phi_deg} deg -> {e}")
        return 0.0
#endfor

def dvcsgen_bh_only(xB, Q2, t_pos, phi_deg, beam_E=10.604, globalfit=True, local=False):
    my_env = os.environ.copy()
    path = REGULAR_DVCSGEN_PATH
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
        line_str = dstot.splitlines()[-1].decode("utf-8")
        tokens = line_str.split()
        val_str = tokens[-1]
        return float(val_str)
    except Exception as e:
        print(f"dvcsgen BH-only error at xB={xB}, Q2={Q2}, t={t_pos}, phi={phi_deg} deg -> {e}")
        return 0.0
#endfor

def dvcsgen_printrad(xB, Q2, t_pos, phi_deg, beam_E=10.604):
    my_env = os.environ.copy()
    path = RADIATIVE_DVCSGEN_PATH

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
        lines = dstot.decode("utf-8", errors="replace").splitlines()
        if len(lines) < 2:
            return (1.0, 0.0)

        penultimate = lines[-2]
        if "Frad_with_error" not in penultimate:
            return (1.0, 0.0)

        tokens = penultimate.split()
        if len(tokens) < 3:
            return (1.0, 0.0)

        factor = float(tokens[1])
        sysval = float(tokens[2])
        return (factor, sysval)
    except Exception as e:
        print(f"dvcsgen printrad error (xB={xB}, Q2={Q2}, t={t_pos}, phi={phi_deg}): {e}")
        return (1.0, 0.0)
#endfor

######################################################
# 2) Parallelizing bin-centering (KM15+VGG) and Frad
######################################################

def _subbin_fbin_task(args):
    """
    Helper function for parallel calls in calculate_fbin.
    Expects args=(xB, Q2, t_pos, beam_E, phi_deg) and returns (km15_val, vgg_val).
    """
    (xB, Q2, t_pos, beam_E, phi_deg) = args
    try:
        km15_val = km15_model(xB, Q2, t_pos, phi_deg, beam_E)
        vgg_val  = dvcsgen_vgg(xB, Q2, t_pos, phi_deg, beam_E, globalfit=False)
        return (km15_val, vgg_val)
    except:
        return (0.0, 0.0)

def _subbin_frad_task(args):
    """
    Helper function for parallel calls in calculate_frad.
    Expects args=(xB, Q2, t_pos, beam_E, phi_deg) => returns single f_sub
    """
    (xB, Q2, t_pos, beam_E, phi_deg) = args
    try:
        f_sub, _subsys = dvcsgen_printrad(xB, Q2, t_pos, phi_deg, beam_E)
        return f_sub
    except:
        return 1.0

######################################################
# 2A) calculate_fbin with concurrency
######################################################
def calculate_fbin(row, prefix, beam_E, n_steps=4):
    """
    Sub-binning with KM15 & VGG from dvcsgen_print code
    to get (km15_fbin, vgg_fbin, final_fbin, fbin_sys_unc).

    We do concurrency to parallelize the sub-bin calls. 
    """

    Mp = 0.938272
    xB_samples   = np.linspace(row['xB_min'], row['xB_max'], n_steps)
    Q2_samples   = np.linspace(row['Q2_min'], row['Q2_max'], n_steps)
    t_pos_samples= np.linspace(row['t_min'],  row['t_max'],  n_steps)
    phi_samples  = np.linspace(row['phi_min'],row['phi_max'],n_steps)

    # We'll gather tasks for each valid sub-bin (w/ cuts) in a list
    tasks = []

    for xB in xB_samples:
        for Q2 in Q2_samples:
            for t_pos in t_pos_samples:
                t_phys = -abs(t_pos)
                try:
                    sqrt_term = np.sqrt(1 + (4*Mp**2*xB**2)/Q2)
                    t_min_val = -Q2*(1 - xB)**2 / (xB*(1+ sqrt_term))
                except:
                    continue

                try:
                    y = Q2/(2*Mp*xB*beam_E)
                    W = np.sqrt(Mp**2 + Q2*(1/xB -1))
                except:
                    continue

                if (t_phys>=t_min_val) and (0.19<y<0.8) and (W>2.0):
                    for phi_deg in phi_samples:
                        tasks.append((xB, Q2, t_pos, beam_E, phi_deg))
                else:
                    continue

    # Parallel execution
    valid_KM15 = []
    valid_VGG  = []

    if len(tasks)==0:
        # fallback if no sub-bins
        pass
    else:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = list(executor.map(_subbin_fbin_task, tasks))
        # results is a list of (km15_val, vgg_val) for each sub-bin
        for (km15_val, vgg_val) in results:
            valid_KM15.append(km15_val)
            valid_VGG.append(vgg_val)

    try:
        centerKM15 = km15_model(row['xB_avg'], row['Q2_avg'], row['t_avg'], row['phi_avg'], beam_E)
        centerVGG = dvcsgen_vgg(row['xB_avg'], row['Q2_avg'], row['t_avg'], row['phi_avg'], beam_E, globalfit=False)
    except:
        return (1.0, 1.0, 1.0, 0.0)

    if not valid_KM15 or not valid_VGG:
        return (1.0, 1.0, 1.0, 0.0)

    avgKM15 = np.mean(valid_KM15)
    avgVGG = np.mean(valid_VGG)
    
    # Calculate standard deviations of sub-bin values
    std_KM15 = np.std(valid_KM15) if len(valid_KM15) > 1 else 0
    std_VGG = np.std(valid_VGG) if len(valid_VGG) > 1 else 0

    # Calculate ratios and their uncertainties
    try:
        km15_fbin = centerKM15 / avgKM15
        # km15_unc = (std_KM15/avgKM15) * km15_fbin  # Relative uncertainty propagation
    except:
        km15_fbin = 1.0
        # km15_unc = 0.0

    try:
        vgg_fbin = centerVGG / avgVGG
        # vgg_unc = (std_VGG/avgVGG) * vgg_fbin
    except:
        vgg_fbin = 1.0
        # vgg_unc = 0.0

    # Calculate model disagreement
    model_disagreement = np.std([km15_fbin, vgg_fbin]) if len(valid_KM15) > 0 and len(valid_VGG) > 0 else 0

    # Combine uncertainties (sub-bin spread + model disagreement)
    total_unc = np.sqrt(
        0 +  # Average sub-bin uncertainty
        model_disagreement**2            # Model-to-model variation
    )

    final_fbin = np.mean([km15_fbin, vgg_fbin])
    
    print("\nUncertainty Breakdown:")
    # print(f"KM15 sub-bin variation: ±{km15_unc:.4f}")
    # print(f"VGG sub-bin variation:  ±{vgg_unc:.4f}")
    print(f"Model disagreement (total systematic):     ±{model_disagreement:.4f}")
    # print(f"Total systematic:       ±{total_unc:.4f}")

    return (km15_fbin, vgg_fbin, final_fbin, total_unc)
#enddef

######################################################
# 1) Radiative sub-binning concurrency tasks
######################################################
def _subbin_frad_task(args):
    """
    Helper function for parallel calls in calculate_frad.
    Expects args=(xB, Q2, t_pos, beam_E, phi_deg) => returns single f_sub
    """
    (xB, Q2, t_pos, beam_E, phi_deg) = args
    try:
        f_sub, _sys_sub = dvcsgen_printrad(xB, Q2, t_pos, phi_deg, beam_E)
        return f_sub
    except:
        return 1.0

######################################################
# 2) calculate_frad with concurrency and clamping (modified)
######################################################
def calculate_frad(row, prefix, beam_E, n_steps=4):
    """
    Modified sub-binning approach returning direct sub-bin average
    Maintains physics validation and error handling
    """
    import math  # For checking isfinite

    print(row)  # debug line

    Mp = 0.938272
    xB_samples = np.linspace(row['xB_min'], row['xB_max'], n_steps)
    Q2_samples = np.linspace(row['Q2_min'], row['Q2_max'], n_steps)
    t_pos_samples = np.linspace(row['t_min'], row['t_max'], n_steps)
    phi_samples = np.linspace(row['phi_min'], row['phi_max'], n_steps)

    tasks = []

    # Physics-validated task generation (keep original validation)
    for xB in xB_samples:
        for Q2 in Q2_samples:
            for t_pos in t_pos_samples:
                t_phys = -abs(t_pos)
                try:
                    sqrt_term = np.sqrt(1 + (4 * Mp**2 * xB**2)/Q2)
                    t_min_val = -Q2*(1 - xB)**2/(xB*(1 + sqrt_term))
                    y = Q2/(2*Mp*xB*beam_E)
                    W = np.sqrt(Mp**2 + Q2*(1/xB -1))
                    
                    if (t_phys >= t_min_val) and (0.19 < y < 0.8) and (W > 2.0):
                        for phi_deg in phi_samples:
                            tasks.append((xB, Q2, t_pos, beam_E, phi_deg))
                except:
                    continue

    subbin_vals = []
    if tasks:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = executor.map(_subbin_frad_task, tasks)
            subbin_vals = list(results)

    # Simplified statistics calculation
    if not subbin_vals:
        return (1.0, 1.0)

    avg_subbin = np.mean(subbin_vals)
    std_subbin = np.std(subbin_vals)

    # Clamping and validation
    if not math.isfinite(avg_subbin) or (avg_subbin < 0) or (avg_subbin > 2):
        return (1.0, 1.0)

    return (avg_subbin, std_subbin)
#enddef

######################################################
# 3) Main code
######################################################
def main():
    n_rows_to_process = parse_args()

    print("Beginning code.")
    input_csv = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/output/unfolding_data.csv"
    print(f"Reading CSV from: {input_csv}")
    df = pd.read_csv(input_csv)
    print(f"Loaded DataFrame with {len(df)} rows and {len(df.columns)} columns.")

    if n_rows_to_process is None:
        n_rows_to_process = len(df)
    else:
        n_rows_to_process = min(n_rows_to_process, len(df))

    fb_cols = ['KM15_Fbin','VGG_Fbin','Fbin','Fbin_sys_uncertainty']
    frad_cols= ['Frad','Frad_sys_uncertainty']

    print("Preparing new columns for Fbin and Frad...")

    for prefix in ['fall','spring']:
        for c in fb_cols:
            df[f"{prefix}_{c}"] = np.nan
        for c in frad_cols:
            df[f"{prefix}_{c}"] = np.nan

    for prefix in ['fall','spring']:
        for c in ['km15','vgg','bh']:
            colname = f"{prefix}_{c}"
            if colname not in df.columns:
                df[colname] = np.nan

    print("Beginning cross-section, Fbin, and Frad calculations...")
    progress_interval = 25

    for i in range(n_rows_to_process):
        if i % progress_interval == 0:
            print(f"  Processing row {i} of {n_rows_to_process} ...")

        row = df.iloc[i]

        # (A) Existing model calcs at bin center
        fall_km15_val = km15_model( row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.604 )
        fall_vgg_val  = dvcsgen_vgg( row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.604, globalfit=False )
        fall_bh_val   = dvcsgen_bh_only( row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.604, globalfit=False )

        spring_km15_val = km15_model( row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.1998 )
        spring_vgg_val  = dvcsgen_vgg( row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.1998, globalfit=False )
        spring_bh_val   = dvcsgen_bh_only( row["xB_avg"], row["Q2_avg"], row["t_avg"], row["phi_avg"], 10.1998, globalfit=False )

        # (B) Fbin for Fall (in parallel)
        fall_KM15_Fbin, fall_VGG_Fbin, fall_Fbin, fall_Fbin_sys = calculate_fbin(row, 'fall', 10.604)
        df.loc[i,"fall_KM15_Fbin"]            = fall_KM15_Fbin
        df.loc[i,"fall_VGG_Fbin"]             = fall_VGG_Fbin
        df.loc[i,"fall_Fbin"]                 = fall_Fbin
        df.loc[i,"fall_Fbin_sys_uncertainty"] = fall_Fbin_sys

        # (C) Frad for Fall => parallel
        fall_Frad_val, fall_Frad_sys = calculate_frad(row, 'fall', 10.604)
        df.loc[i,"fall_Frad"] = fall_Frad_val
        df.loc[i,"fall_Frad_sys_uncertainty"] = fall_Frad_sys

        # (D) Fbin for Spring
        spring_KM15_Fbin, spring_VGG_Fbin, spring_Fbin, spring_Fbin_sys = calculate_fbin(row, 'spring', 10.1998)
        df.loc[i,"spring_KM15_Fbin"]            = spring_KM15_Fbin
        df.loc[i,"spring_VGG_Fbin"]             = spring_VGG_Fbin
        df.loc[i,"spring_Fbin"]                 = spring_Fbin
        df.loc[i,"spring_Fbin_sys_uncertainty"] = spring_Fbin_sys

        # (E) Frad for Spring
        spring_Frad_val, spring_Frad_sys = calculate_frad(row, 'spring', 10.1998)
        df.loc[i,"spring_Frad"]                 = spring_Frad_val
        df.loc[i,"spring_Frad_sys_uncertainty"] = spring_Frad_sys

        # (F) Store bin-center model outputs
        df.loc[i,"fall_km15"]   = fall_km15_val
        df.loc[i,"fall_vgg"]    = fall_vgg_val
        df.loc[i,"fall_bh"]     = fall_bh_val
        df.loc[i,"spring_km15"] = spring_km15_val
        df.loc[i,"spring_vgg"]  = spring_vgg_val
        df.loc[i,"spring_bh"]   = spring_bh_val

        # (G) Apply final (Fbin * Frad) to cross section & stat unc
        for pre in ['fall','spring']:
            fbin_factor = df.loc[i,f"{pre}_Fbin"]
            frad_factor= df.loc[i,f"{pre}_Frad"]
            total_factor = fbin_factor*frad_factor
            df.loc[i,f"{pre}_cross_section"] *= total_factor
            df.loc[i,f"{pre}_cross_section_stat_uncertainty"] *= total_factor

    print("Reordering columns...")
    all_cols = list(df.columns)

    def get_ordered_columns(prefix):
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
    common_columns = [c for c in all_cols if not c.startswith(("fall_","spring_"))]

    fall_cols = get_ordered_columns("fall")
    spring_cols = get_ordered_columns("spring")

    for col in common_columns:
        new_column_order.append(col)
        if col in all_cols:
            all_cols.remove(col)

    for col in fall_cols:
        if col in all_cols:
            new_column_order.append(col)
            all_cols.remove(col)

    for col in spring_cols:
        if col in all_cols:
            new_column_order.append(col)
            all_cols.remove(col)

    new_column_order += all_cols
    df = df[new_column_order]

    outpath = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/output/unfolding_data_with_models.csv"
    df.to_csv(outpath, index=False)
    print(f"Done! Updated file saved to:\n{outpath}")


if __name__ == "__main__":
    main()