// Third panel: x histograms scaled by the fit constant with propagated errors
    c1->cd(1);
    gPad->SetLeftMargin(0.15);
    TH1D *h_x_nh3 = 
        new TH1D("h_x_nh3", "P_{T} Distribution; P_{T} (GeV); Counts", 100, 0.05, 1.0);
    TH1D *h_x_carbon = 
        new TH1D("h_x_carbon", "P_{T} Distribution; P_{T} (GeV); Counts", 100, 0.05, 1.0);
    tree_nh3->Draw("x>>h_x_nh3");
    tree_carbon->Draw("x>>h_x_carbon");
    TH1D *h_x_carbon_scaled = (TH1D*)h_x_carbon->Clone("h_x_carbon_scaled");
    h_x_carbon_scaled->SetTitle("x_B Distribution; x_B; Counts (Scaled)");

    for (int i = 1; i <= h_x_carbon->GetNbinsX(); ++i) {
        double bin_content = h_x_carbon->GetBinContent(i);
        double bin_error = h_x_carbon->GetBinError(i);

        double new_content = bin_content * scale_factor;
        double new_error = 
            new_content * std::sqrt((bin_error / bin_content) * (bin_error / bin_content) + 
            (scale_error / scale_factor) * (scale_error / scale_factor));
        h_x_carbon_scaled->SetBinContent(i, new_content);
        h_x_carbon_scaled->SetBinError(i, new_error);
    }

    TGraphErrors *gr_dilution_x = new TGraphErrors();
    for (int i = 1; i <= h_x_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_x_nh3->GetBinContent(i);
        double nh3_error = h_x_nh3->GetBinError(i);
        double c_counts = h_x_carbon_scaled->GetBinContent(i);
        double c_error = h_x_carbon_scaled->GetBinError(i);

        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;

            // Propagate the error
            double dilution_error = std::sqrt(
                std::pow((c_counts / (nh3_counts * nh3_counts)) * nh3_error, 2) +
                std::pow(1.0 / nh3_counts * c_error, 2));

            gr_dilution->SetPoint(i - 1, h_x_nh3->GetBinCenter(i), dilution);
            gr_dilution->SetPointError(i - 1, 0, dilution_error);
        }
    }
    gr_dilution_x->SetTitle("; x_B; D_f = (NH3 - s*C) / NH3");
    gr_dilution_x->SetMarkerStyle(20);
    gr_dilution_x->Draw("AP");
    gr_dilution_x->GetXaxis()->SetRangeUser(0, 1);
    gr_dilution_x->GetYaxis()->SetRangeUser(0.05, 0.15);

    // Fit to a third-degree polynomial
    TF1 *fit_poly_x = new TF1("fit_poly", "[0] + [1]*x + [2]*x^2 + [3]*x^3", 0.05, 1.0);
    gr_dilution_x->Fit(fit_poly, "R");
    fit_poly_x->SetLineColor(kRed);
    fit_poly_x->Draw("SAME");

    // Retrieve fit parameters and their errors
    p0_x = fit_poly->GetParameter(0);
    p0_err_x = fit_poly->GetParError(0);
    p1_x = fit_poly->GetParameter(1);
    p1_err_x = fit_poly->GetParError(1);
    p2_x = fit_poly->GetParameter(2);
    p2_err_x = fit_poly->GetParError(2);
    p3_x = fit_poly->GetParameter(3);
    p3_err_x = fit_poly->GetParError(3);

    // Retrieve chi2 and NDF
    chi2 = fit_poly->GetChisquare();
    ndf = fit_poly->GetNDF();
    chi2_ndf = chi2 / ndf;

    // Add fit parameters box
    TPaveText *x = new TPaveText(0.5, 0.7, 0.9, 0.9, "brNDC");
    x->SetBorderSize(1);
    x->SetFillStyle(1001); // Solid fill style
    x->SetFillColor(kWhite); // White background
    x->AddText(Form("p0 = %.3f +/- %.3f", p0, p0_err));
    x->AddText(Form("p1 = %.3f +/- %.3f", p1, p1_err));
    x->AddText(Form("p2 = %.3f +/- %.3f", p2, p2_err));
    x->AddText(Form("p3 = %.3f +/- %.3f", p3, p3_err));
    x->Draw();

    // Add chi2/ndf in the top left
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.15, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2, ndf, chi2_ndf));