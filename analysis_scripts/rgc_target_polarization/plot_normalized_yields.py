def plot_normalized_yields(trees, xB_bins):
    """
    For each target = NH3, C, CH2, He, ET, plot the x_B distribution
    for RGC_Su22, RGC_Fa22 and RGC_Sp23 normalized by total accumulated charge.
    Arranged on a 2x3 canvas; the 6th subplot is left blank.
    Saves to output/normalized_yields.pdf.
    """
    # absolute accumulated charges (in same units as counts denominator)
    CHARGE = {
        "RGC_Su22": {"NH3": 3686969.636627, "C": 363715.413199, "CH2": 189723.230200, "He": 382585.145820, "ET": 470011.187943},
        "RGC_Fa22": {"NH3": 5509572.178076, "C":2006703.362768, "CH2":1833777.957336, "He":378165.719362, "ET": 57332.748981},
        "RGC_Sp23": {"NH3":1948649.296492, "C":552805.225420, "CH2":436816.389755, "He":279407.815704, "ET":171738.938831},
    }

    periods = ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]
    targets = ["NH3", "C", "CH2", "He", "ET"]

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, target in enumerate(targets):
        ax = axes[i]
        for period in periods:
            tree = trees[period][target]
            # pull x array
            x = tree["x"].array(library="np")
            # histogram counts & edges
            counts, edges = np.histogram(x, bins=xB_bins)
            centers = 0.5*(edges[:-1] + edges[1:])
            # normalize by total charge
            norm_counts = counts / CHARGE[period][target]
            ax.step(centers, norm_counts, where="mid", label=period)
        ax.set_xlabel("$x_{B}$")
        ax.set_ylabel("counts/nC")
        ax.set_title(target)
        ax.legend(fontsize="small")
        ax.grid(True)

    # turn off the unused 6th subplot
    axes[-1].axis("off")

    plt.tight_layout()
    os.makedirs("output", exist_ok=True)
    out_path = "output/normalized_yields.pdf"
    fig.savefig(out_path)
    plt.close(fig)
    print(f"[Plot] Saved normalized yields to {out_path}")