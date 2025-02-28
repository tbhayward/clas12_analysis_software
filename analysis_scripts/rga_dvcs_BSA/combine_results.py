def combine_results(output_dir):
    """
    Optional final JSON combiner
    """
    import json

    combined = {}
    all_periods = [
        "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb",
        "eppi0_Fa18_inb", "eppi0_Fa18_out", "eppi0_Sp19_inb"
    ]
    topologies = ["FD_FD", "CD_FD", "CD_FT"]

    for period_code in all_periods:
        for topo in topologies:
            fname = f"cuts_{period_code}_{topo}_final.json"
            fpath = os.path.join(output_dir, fname)
            if not os.path.exists(fpath):
                continue
            #endif
            key_name = f"{period_code}_{topo}"
            try:
                with open(fpath, "r") as f:
                    combined[key_name] = json.load(f)
                #endwith
            except FileNotFoundError:
                print(f"⚠️ Missing file {fpath}")
            #endtry
        #endfor
    #endfor

    combined_path = os.path.join(output_dir, "combined_cuts.json")
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2)
    #endif
    print(f"✅ Wrote combined JSON => {combined_path}")
#enddef