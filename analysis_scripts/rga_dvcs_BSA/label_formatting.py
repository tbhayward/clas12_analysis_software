# label_formatting.py

import ROOT

def configure_global_style():
    """
    Applies a global ROOT style for consistent plotting.
    Called once at the beginning in main.py.
    """
    style = ROOT.gStyle
    style.SetOptStat(0)
    style.SetTitleSize(0.045, "XY")
    style.SetLabelSize(0.04, "XY")
    style.SetPadLeftMargin(0.12)
    style.SetPadRightMargin(0.08)
    style.SetPadTopMargin(0.08)
    style.SetPadBottomMargin(0.12)
    style.SetLegendBorderSize(1)
    style.SetLegendTextSize(0.035)
    style.SetLegendFillColor(0)
    style.SetGridColor(ROOT.kGray+1)
    style.SetGridStyle(2)
    style.SetGridWidth(1)
#enddef

def format_label_name(variable, analysis_type):
    """
    Returns axis label text for a given variable & analysis type.
    If 'analysis_type' == 'eppi0', we might replace gamma with pi0 in the label.
    """
    labels = {
        "open_angle_ep2": "#theta_{e'#gamma} (deg)",
        "Delta_phi": "#Delta#phi (rad)",
        "theta_gamma_gamma": "#theta_{#gamma#gamma} (deg)",
        "theta_pi0_pi0": "#theta_{#pi^{0}#pi^{0}} (deg)",
        "Emiss2": "Miss E^{2} (GeV^{2})",
        "Mx2": "Miss M^{2} (GeV^{2})",
        "Mx2_1": "Miss M_{e'p'}^{2} (GeV^{2})",
        "Mx2_2": "Miss M_{e'#gamma}^{2} (GeV^{2})",
        "pTmiss": "p_{T}^{miss} (GeV)",
        "xF": "x_{F}",
        # We can also do an open_angle_ep2 if eppi0 vs dvcs, etc.
    }

    if analysis_type == "eppi0":
        # eppi0 might need to replace some gamma references with pi0.
        # For instance, we might rename Mx2_2 to Miss M_{e'#pi^{0}}^{2}, etc.
        # But below is just an example:
        return labels.get(variable, variable).replace("#gamma", "#pi^{0}")
    else:
        return labels.get(variable, variable)
    #endif
#enddef