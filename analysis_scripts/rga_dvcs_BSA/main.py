# main.py

import os

# Import only the pieces needed to run your full analysis:
from exclusivity import process_period, combine_results
# If you want to configure global style, import it from label_formatting or wherever it is:
from label_formatting import configure_global_style

def main():
    """
    Master workflow for the analysis.
    - Configure style (if you like).
    - Create output directory.
    - Loop over run periods.
    - Call process_period(...) for each.
    - Combine final JSON outputs.
    """
    configure_global_style()  # Or call it once, if needed

    output_dir = "exclusivity_plots"
    os.makedirs(output_dir, exist_ok=True)

    print("🚀 Starting analysis with sequential processing\n")

    # Example of processing multiple run periods
    for period in ["Fa18_inb", "Fa18_out", "Sp19_inb"]:
        process_period(period, output_dir)
    #endfor

    print("🧩 Combining results...")
    combine_results(output_dir)
    print("\n🎉 Analysis complete!")
#enddef


if __name__ == "__main__":
    main()
#endif