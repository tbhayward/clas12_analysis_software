# kin_cuts.py

def apply_kinematic_cuts(t_value, open_angle_ep2_value, theta_neutral_neutral_value,
                         Emiss2_value, Mx2_value, Mx2_1_value, Mx2_2_value,
                         pTmiss_value, xF_value, analysis_type, data_type,
                         run_period, topology):
    """
    A placeholder for real kinematic cuts, but we do enforce open_angle_ep2 > 5
    from the very beginning. This ensures we never fill histograms for events with
    open_angle_ep2 <= 5, even in the first pass.

    If you want to add more universal cuts, do so here.
    """
    # 1) universal angle cut
    if open_angle_ep2_value <= 5.0:
        return False
    #endif

    # 2) universal angle cut
    if (-t_value) > 1.0:
        return False
    #endif

    # 3) universal missing pt cut
    if pTmiss_value > 0.20:
        return False
    #endif

    # If all checks pass:
    print("We're going to return true now!");
    return True
#enddef


def passes_3sigma_cuts(event, is_mc, cuts_dict):
    """
    Enforce a mu +/- 3*sigma cut for each variable in 'cuts_dict'.
    'cuts_dict' looks like:
      {
        "data": { var_name: {"mean": X, "std": Y}, ... },
        "mc":   { var_name: {"mean": X, "std": Y}, ... }
      }
    If event.var is outside mu +/- 3*sigma for any variable, we return False.
    """
    category = "mc" if is_mc else "data"

    # Iterate over all variables we have in this dictionary:
    for var, stats in cuts_dict.get(category, {}).items():
        mean = stats["mean"]
        stdv = stats["std"]
        val = getattr(event, var, None)
        if val is None:
            # If variable doesn't exist in the event, skip or keep? We'll skip.
            continue
        #endif

        # The actual 2.5sigma check:
        if val < (mean - 2.5*stdv) or val > (mean + 2.5*stdv):
            return False
        #endif
    #endfor

    return True
#enddef