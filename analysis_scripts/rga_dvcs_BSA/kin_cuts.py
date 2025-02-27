# kin_cuts.py

def apply_kinematic_cuts(t_value, open_angle_ep2_value, theta_neutral_neutral_value,
                         Emiss2_value, Mx2_value, Mx2_1_value, Mx2_2_value, pTmiss_value,
                         xF_value, analysis_type, data_type, run_period, topology):
    """
    Placeholder (or real) kinematic cuts logic.
    Currently always returns True.
    """
    return True
#enddef

def passes_3sigma_cuts(event, is_mc, cuts_dict):
    """
    Enforce a mu +/- 3*sigma cut for each variable in the 'cuts_dict'.
    'cuts_dict' is loaded from JSON elsewhere (but used here).
    """
    category = "mc" if is_mc else "data"
    for var, stats in cuts_dict[category].items():
        mean = stats["mean"]
        stdv = stats["std"]
        val = getattr(event, var, None)
        if val is None:
            # If variable not found in event, skip or handle gracefully
            continue
        #endif
        if val < (mean - 3.0*stdv) or val > (mean + 3.0*stdv):
            return False
        #endif
    #endfor
    return True
#enddef