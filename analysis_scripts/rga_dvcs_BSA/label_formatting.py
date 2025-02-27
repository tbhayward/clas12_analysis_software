# label_formatting.py

import ROOT
# Configure ROOT to run in batch mode
ROOT.gROOT.SetBatch(True)

def configure_global_style():
    """
    Applies a global ROOT style for consistent plotting.
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
    """
    labels = {
        "open_angle_ep2": "#theta_{e'#gamma} (deg)",
        "Mx2_2": "Miss M_{e'#gamma}^{2} (GeV^{2})",
        "theta_gamma_gamma": "#theta_{#gamma#gamma} (deg)",
        "theta_pi0_pi0": "#theta_{#pi^{0}#pi^{0}} (deg)",
        "xF": "x_{F}",
        "Emiss2": "Miss E^{2} (GeV^{2})",
        "Mx2": "Miss M^{2} (GeV^{2})",
        "Mx2_1": "Miss M_{e'p'}^{2} (GeV^{2})",
        "pTmiss": "p_{T}^{miss} (GeV)"
    }
    if analysis_type == "eppi0":
        return labels.get(variable, variable).replace("#gamma", "#pi^{0}")
    return labels.get(variable, variable)
#enddef