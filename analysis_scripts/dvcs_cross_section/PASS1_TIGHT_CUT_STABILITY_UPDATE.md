# Pass-1 tight-cut stability prescription

The exclusivity- and fiducial-cut systematic calculations now use the pass-1
stability prescription in each bin.

For nominal, loose and tight cross sections, define

    d_loose = sigma_loose - sigma_nominal
    d_tight = sigma_tight - sigma_nominal

If

    |d_tight| / |sigma_nominal| > 0.50,

the tight result is treated as unstable because of insufficient statistics and
only the loose variation is used:

    s_cut = |d_loose|.

Otherwise,

    s_cut = 0.5 * (|d_loose| + |d_tight|).

The stored systematic remains an absolute uncertainty in nb/GeV^4.

Barlow filtering remains active. The `raw` column applies the prescription to
both unfiltered deviations. The final column first sets each deviation that
fails B >= 1 to zero, then applies the same pass-1 prescription. The decision
to invoke the loose-only branch is based on the unfiltered tight relative
difference, exactly matching the stated 50% stability criterion.

The diagnostic CSV now includes, separately for exclusivity and fiducial cuts:

- tight relative difference
- whether the loose-only branch was used
- loose and tight Barlow values
- loose and tight retained flags
- raw and final absolute systematics

The automatic runner exposes these settings in `main.cpp`:

    cut_variation_opts.use_pass1_tight_instability_rule = true;
    cut_variation_opts.tight_relative_difference_threshold = 0.50;
