# main.py

import os

# We import the multi-stage process_period function and the combine_results function
from exclusivity import process_period_multi_stage, combine_results
# For global ROOT styling (if you desire):
from label_formatting import configure_global_style

def main():
    """
    The 'main' entry point.
    - Applies global style (optional).
    - Creates output directory.
    - Loops over run periods (DVCS or eppi0).
    - Calls process_period_multi_stage(...) for each.
    - Combines final JSON outputs if you like.
    """
    # Optional: configure global ROOT style once at startup
    configure_global_style()
    
    output_dir = "exclusivity"
    os.makedirs(output_dir, exist_ok=True)

    print("🚀 Starting multi-stage analysis with sequential processing\n")

    # EXAMPLE: Suppose you want to process DVCS for Fa18_inb, Fa18_out, Sp19_inb,
    # and also eppi0 for some set of periods. You can specify them:
    periods_to_run = [
        ("DVCS_Fa18_inb",  "dvcs"),
        ("DVCS_Fa18_out",  "dvcs"),
        ("DVCS_Sp19_inb",  "dvcs"),
        ("eppi0_Fa18_inb", "eppi0"),
        ("eppi0_Fa18_out", "eppi0"),
        ("eppi0_Sp19_inb", "eppi0"),
    ]

    for (period, analysis_type) in periods_to_run:
        process_period_multi_stage(period, output_dir, analysis_type)
    #endfor

    print("🧩 Combining results (JSON files from each topology & stage)...")
    combine_results(output_dir)
    print("\n🎉 Analysis complete!")
#enddef

if __name__ == "__main__":
    main()
#endif