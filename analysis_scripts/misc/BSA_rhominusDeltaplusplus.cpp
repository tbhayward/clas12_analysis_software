// Created 9/25/23
// get BSA of rho- Delta++ 

#include <TCanvas.h>
#include <TH1F.h>
#include <TObjArray.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <iostream>
#include <string>
#include <map>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

// Function to set the style of a canvas and its pads
void setCanvasStyle(TCanvas &canvas, int nCols, int nRows) {
    canvas.Divide(nCols, nRows);
    for (int i = 1; i <= nCols * nRows; ++i) {
        canvas.cd(i);
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.12);
    }
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> calculateAndPlotALU(
    TTreeReader &dataReader, const char* branchName, double min_val, double max_val,
    int num_kinematic_bins, int fit) {

    std::vector<double> ALU_values;
    std::vector<double> ALU_errors;

    // Create a 2D array to hold N+ and N- for each dynamic bin and phi bin
    std::vector<std::vector<double>> N_pos(num_kinematic_bins, std::vector<double>(12, 0));  
    // "num_kinematic_bins" dynamic bins, 12 phi bins
    std::vector<std::vector<double>> N_neg(num_kinematic_bins, std::vector<double>(12, 0));

    std::vector<double> sum_beam_pol(num_kinematic_bins, 0.0);
    std::vector<int> count_beam_pol(num_kinematic_bins, 0);
    std::vector<double> sum_W_over_A(num_kinematic_bins, 0.0);
    std::vector<int> count_W_over_A(num_kinematic_bins, 0);
    // Declare additional vectors to hold the sum and count of each dynamic bin.
    std::vector<double> sum_branch_var(num_kinematic_bins, 0.0);
    std::vector<int> count_branch_var(num_kinematic_bins, 0);


    // Declare reader locations
    TTreeReaderValue<int> runnum(dataReader, "runnum");
    TTreeReaderValue<int> helicity(dataReader, "helicity");
    TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
    TTreeReaderValue<double> e_p(dataReader, "e_p");
    TTreeReaderValue<double> e_theta(dataReader, "e_theta");
    TTreeReaderValue<double> e_phi(dataReader, "e_phi");
    TTreeReaderValue<double> p1_p(dataReader, "p1_p");
    TTreeReaderValue<double> p1_theta(dataReader, "p1_theta");
    TTreeReaderValue<double> p1_phi(dataReader, "p1_phi");
    TTreeReaderValue<double> p2_p(dataReader, "p2_p");
    TTreeReaderValue<double> p2_theta(dataReader, "p2_theta");
    TTreeReaderValue<double> p2_phi(dataReader, "p2_phi");
    TTreeReaderValue<double> p3_p(dataReader, "p3_p");
    TTreeReaderValue<double> p3_theta(dataReader, "p3_theta");
    TTreeReaderValue<double> p3_phi(dataReader, "p3_phi");
    TTreeReaderValue<double> phi1(dataReader, "phi1");
    TTreeReaderValue<double> phi2(dataReader, "phi2");
    TTreeReaderValue<double> phi3(dataReader, "phi3");
    TTreeReaderValue<double> phi12(dataReader, "phi12");
    std::unique_ptr<TTreeReaderValue<double>> phi;

    if (fit == 0) {
        phi = std::make_unique<TTreeReaderValue<double>>(dataReader, "phi13");
    } else {
        phi = std::make_unique<TTreeReaderValue<double>>(dataReader, "phi23");
    }

    TTreeReaderValue<double> DepW(dataReader, "DepW");
    TTreeReaderValue<double> DepA(dataReader, "DepA");
    TTreeReaderValue<double> x(dataReader, "x");
    TTreeReaderValue<double> z1(dataReader, "z1");
    TTreeReaderValue<double> z13(dataReader, "z13");
    TTreeReaderValue<double> Mx(dataReader, "Mx");
    TTreeReaderValue<double> Mx12(dataReader, "Mx12");
    TTreeReaderValue<double> Mx13(dataReader, "Mx13");
    TTreeReaderValue<double> Mx23(dataReader, "Mx23");
    TTreeReaderValue<double> Mh12(dataReader, "Mh12");
    TTreeReaderValue<double> branch_var(dataReader, branchName);
    TTreeReaderValue<double> Mh23(dataReader, "Mh23");
    TTreeReaderValue<double> xF13(dataReader, "xF13");
    TTreeReaderValue<double> Delta_phi13(dataReader, "Delta_phi13");

    // Declare new variables to store missing particle information
    float px_p, px_theta, px_phi, Mx1x, Mx2x, Mx3x, Mh1x, Mh2x, Mh3x;

    // Define initial state 4-momentum (10.1998 GeV electron beam and stationary proton)
    TLorentzVector p_initial(0, 0, 10.1998, 10.1998 + 0.938); // (px, py, pz, E)

    int counter = 0;
    while (dataReader.Next()) {
        counter++;
        // if (counter > 1000000) { break; }
        if (*Mx > 0.3) continue;

        // Create 4-momentum vectors for final state particles
        TLorentzVector p_e, p1, p2, p3;
        p_e.SetXYZM(*e_p*sin(*e_theta)*cos(*e_phi), *e_p*sin(*e_theta)*sin(*e_phi), 
            *e_p*cos(*e_theta), 0.511e-3);
        p1.SetXYZM(*p1_p*sin(*p1_theta)*cos(*p1_phi), *p1_p*sin(*p1_theta)*sin(*p1_phi), 
            *p1_p*cos(*p1_theta), 0.139570);
        p2.SetXYZM(*p2_p*sin(*p2_theta)*cos(*p2_phi), *p2_p*sin(*p2_theta)*sin(*p2_phi), 
            *p2_p*cos(*p2_theta), 0.139570);
        p3.SetXYZM(*p3_p*sin(*p3_theta)*cos(*p3_phi), *p3_p*sin(*p3_theta)*sin(*p3_phi), 
            *p3_p*cos(*p3_theta), 0.938272);
        // Calculate 4-momentum of missing particle
        TLorentzVector p_x = p_initial - (p_e + p1 + p2 + p3);
        // Populate missing particle variables
        px_p = p_x.P();
        px_theta = p_x.Theta();
        px_phi = p_x.Phi();
        // Calculate missing mass variables
        Mx1x = (p_initial - (p_e + p1)).M();
        Mx2x = (p_initial - (p_e + p2)).M();
        Mx3x = (p_initial - (p_e + p3)).M();
        // Calculate invariant mass variables
        Mh1x = (p1 + p_x).M();
        Mh2x = (p2 + p_x).M();
        Mh3x = (p3 + p_x).M();

        if (Mh2x < 0.775 - 0.13 || Mh2x > 0.775 + 0.13) continue;

        if(*branch_var < min_val || *branch_var > max_val) continue;  
        // Skip entries out of range
        if(**phi < 0 || **phi > 2 * TMath::Pi()) continue;  
        // Skip entries out of range

        int dyn_bin = int((*branch_var - min_val) / ((max_val - min_val) / num_kinematic_bins));
        int phi_bin = int(**phi / (2 * TMath::Pi() / 12));

        if(dyn_bin < 0 || dyn_bin >= num_kinematic_bins) continue;  // Skip invalid indices
        if(phi_bin < 0 || phi_bin >= 12) continue;  // Skip invalid indices

        if (*helicity > 0) {
            N_pos[dyn_bin][phi_bin]++;
        } else if (*helicity < 0) {
            N_neg[dyn_bin][phi_bin]++;
        }
        if (*helicity != 0) {  // assuming non-zero helicity implies valid polarization
            sum_beam_pol[dyn_bin] += *beam_pol;
            count_beam_pol[dyn_bin]++;

            double W_over_A = (*DepA != 0) ? *DepW / *DepA : 0;
            sum_W_over_A[dyn_bin] += W_over_A;
            count_W_over_A[dyn_bin]++;

            sum_branch_var[dyn_bin] += *branch_var;
            count_branch_var[dyn_bin]++;
        }
    }

    for (int dyn_bin = 0; dyn_bin < num_kinematic_bins; ++dyn_bin) {
        TF1 fitFunc("fitFunc", "[0]*sin(x)", 0, 2 * TMath::Pi());
        TGraphErrors fitGraph;
        for (int phi_bin = 0; phi_bin < 12; ++phi_bin) {
            double phi_val = phi_bin * (2 * TMath::Pi() / 12);
            double mean_beam_pol = (count_beam_pol[dyn_bin] != 0) ? sum_beam_pol[dyn_bin] 
                / count_beam_pol[dyn_bin] : 1.0;
            double ALU = (1 / mean_beam_pol) * 
                (N_pos[dyn_bin][phi_bin] - N_neg[dyn_bin][phi_bin]) / 
                (N_pos[dyn_bin][phi_bin] + N_neg[dyn_bin][phi_bin]);
            double ALU_error = (2 / mean_beam_pol) * 
                TMath::Sqrt((N_pos[dyn_bin][phi_bin] * N_neg[dyn_bin][phi_bin]) / 
                TMath::Power(N_pos[dyn_bin][phi_bin] + N_neg[dyn_bin][phi_bin], 3));

            fitGraph.SetPoint(phi_bin, phi_val, ALU);
            fitGraph.SetPointError(phi_bin, 0, ALU_error);
        }
        fitGraph.Fit(&fitFunc, "Q");
        double A = fitFunc.GetParameter(0);
        double A_error = fitFunc.GetParError(0);
        double mean_W_over_A = (count_W_over_A[dyn_bin] != 0) ? 
            sum_W_over_A[dyn_bin] / count_W_over_A[dyn_bin] : 1.0;
        ALU_values.push_back(A / mean_W_over_A);
        ALU_errors.push_back(A_error / mean_W_over_A);
    }   

    // Declare a new vector to hold the average bin values.
    std::vector<double> average_bin_values(num_kinematic_bins, 0.0);
    // Calculate the average bin values.
    for (int i = 0; i < num_kinematic_bins; ++i) {
        if (count_branch_var[i] != 0) {
            average_bin_values[i] = sum_branch_var[i] / count_branch_var[i];
        } else {
            average_bin_values[i] = 0.0;
        }
    }
    // Return all the calculated values and errors, including the average bin values.
    return std::make_tuple(ALU_values, ALU_errors, average_bin_values);
}

void createBSAPlot(TTreeReader &dataReader, const char* outDir) {
    // Declare derived variables to read from the tree
    double z2x, t13, t2x;
    // Declare reader locations
    TTreeReaderValue<int> helicity(dataReader, "helicity");
    TTreeReaderValue<double> beam_pol(dataReader, "beam_pol");
    TTreeReaderValue<double> e_p(dataReader, "e_p");
    TTreeReaderValue<double> e_theta(dataReader, "e_theta");
    TTreeReaderValue<double> e_phi(dataReader, "e_phi");
    TTreeReaderValue<double> p1_p(dataReader, "p1_p");
    TTreeReaderValue<double> p1_theta(dataReader, "p1_theta");
    TTreeReaderValue<double> p1_phi(dataReader, "p1_phi");
    TTreeReaderValue<double> p2_p(dataReader, "p2_p");
    TTreeReaderValue<double> p2_theta(dataReader, "p2_theta");
    TTreeReaderValue<double> p2_phi(dataReader, "p2_phi");
    TTreeReaderValue<double> p3_p(dataReader, "p3_p");
    TTreeReaderValue<double> p3_theta(dataReader, "p3_theta");
    TTreeReaderValue<double> p3_phi(dataReader, "p3_phi");
    TTreeReaderValue<double> phi1(dataReader, "phi1");
    TTreeReaderValue<double> phi2(dataReader, "phi2");
    TTreeReaderValue<double> phi3(dataReader, "phi3");
    TTreeReaderValue<double> phi12(dataReader, "phi12");
    TTreeReaderValue<double> phi13(dataReader, "phi13");
    TTreeReaderValue<double> phi23(dataReader, "phi23");
    TTreeReaderValue<double> DepW(dataReader, "DepW");
    TTreeReaderValue<double> DepA(dataReader, "DepA");
    TTreeReaderValue<double> x(dataReader, "x");
    TTreeReaderValue<double> z1(dataReader, "z1");
    TTreeReaderValue<double> z13(dataReader, "z13");
    TTreeReaderValue<double> Mx(dataReader, "Mx");
    TTreeReaderValue<double> Mx12(dataReader, "Mx12");
    TTreeReaderValue<double> Mx13(dataReader, "Mx13");
    TTreeReaderValue<double> Mx23(dataReader, "Mx23");
    TTreeReaderValue<double> Mh12(dataReader, "Mh12");
    TTreeReaderValue<double> Mh13(dataReader, "Mh13");
    TTreeReaderValue<double> Mh23(dataReader, "Mh23");
    TTreeReaderValue<double> xF13(dataReader, "xF13");
    TTreeReaderValue<double> Delta_phi13(dataReader, "Delta_phi13");

    // Declare new variables to store missing particle information
    float px_p, px_theta, px_phi, Mx1x, Mx2x, Mx3x, Mh1x, Mh2x, Mh3x;

    // Define initial state 4-momentum (10.1998 GeV electron beam and stationary proton)
    TLorentzVector p_initial(0, 0, 10.1998, 10.1998 + 0.938); // (px, py, pz, E)

    // Create canvas and set its style
    TCanvas canvas("Asymmetry", "Canvas", 1600, 1000);
    setCanvasStyle(canvas, 1, 2);

    // Declare a temporary histogram to get statistics
    TH1F tempHist("bin hist", "", 1000, 1, 2.2);

    int counter = 0;
    while (dataReader.Next()) {
        counter++;
        if (counter > 100000) { break; }
        if (*Mx < 0 || *Mx12 < 0 || *Mx13 < 0 || *Mx23 < 0) { continue; }
        if (*Mx > 0.30) { continue; }
        // Create 4-momentum vectors for final state particles
        TLorentzVector p_e, p1, p2, p3;
        p_e.SetXYZM(*e_p*sin(*e_theta)*cos(*e_phi), *e_p*sin(*e_theta)*sin(*e_phi), 
            *e_p*cos(*e_theta), 0.511e-3);
        p1.SetXYZM(*p1_p*sin(*p1_theta)*cos(*p1_phi), *p1_p*sin(*p1_theta)*sin(*p1_phi), 
            *p1_p*cos(*p1_theta), 0.139570);
        p2.SetXYZM(*p2_p*sin(*p2_theta)*cos(*p2_phi), *p2_p*sin(*p2_theta)*sin(*p2_phi), 
            *p2_p*cos(*p2_theta), 0.139570);
        p3.SetXYZM(*p3_p*sin(*p3_theta)*cos(*p3_phi), *p3_p*sin(*p3_theta)*sin(*p3_phi), 
            *p3_p*cos(*p3_theta), 0.938272);
        // Calculate 4-momentum of missing particle
        TLorentzVector p_x = p_initial - (p_e + p1 + p2 + p3);
        // Populate missing particle variables
        px_p = p_x.P();
        px_theta = p_x.Theta();
        px_phi = p_x.Phi();
        // Calculate missing mass variables
        Mx1x = (p_initial - (p_e + p1)).M();
        Mx2x = (p_initial - (p_e + p2)).M();
        Mx3x = (p_initial - (p_e + p3)).M();
        // Calculate invariant mass variables
        Mh1x = (p1 + p_x).M();
        Mh2x = (p2 + p_x).M();
        Mh3x = (p3 + p_x).M();

        if (Mh2x > (0.775 - 0.13) && Mh2x < (0.775 + 0.13)) { tempHist.Fill(*Mh13); }
    }
    dataReader.Restart();  // Reset the TTreeReader at the end of the function

    // Find the quantile edges
    int nQuantiles = 9;
    double quantiles[nQuantiles];
    double sum = tempHist.GetEntries();
    for (int i = 1; i <= nQuantiles; ++i) {
        quantiles[i-1] = i * (sum / nQuantiles);
    }
    double edges[nQuantiles + 1];
    tempHist.GetQuantiles(nQuantiles, edges, quantiles);
    // Get min and max values for the branch
    double min_val = tempHist.GetXaxis()->GetXmin();
    double max_val = tempHist.GetXaxis()->GetXmax();

    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> resultMh13;
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> resultMh23;
    int num_kinematic_bins = 20;
    resultMh13 = calculateAndPlotALU(dataReader, "Mh13", min_val, max_val, num_kinematic_bins, 0);
    dataReader.Restart();  // Reset the TTreeReader at the end of the function
    resultMh23 = calculateAndPlotALU(dataReader, "Mh23", min_val, max_val, num_kinematic_bins, 1);

    // Extract the average_bin_values
    std::vector<double> average_bin_values1 = std::get<2>(resultMh13);
    std::vector<double> average_bin_values2 = std::get<2>(resultMh23);
    TGraphErrors aluGraph1(num_kinematic_bins);
    TGraphErrors aluGraph2(num_kinematic_bins);

    double offset = 0.1 * ((max_val - min_val) / num_kinematic_bins);  // 10% of bin width
    // Populate aluGraph using result
    for (int dyn_bin = 0; dyn_bin < num_kinematic_bins; ++dyn_bin) {
        if (std::get<0>(resultMh13)[dyn_bin] == 0 ) { continue; }
        // Use average_bin_values instead of bin_center
        aluGraph1.SetPoint(dyn_bin, average_bin_values1[dyn_bin], 
            std::get<0>(resultMh13)[dyn_bin]);
        aluGraph1.SetPointError(dyn_bin, 0, std::get<1>(resultMh13)[dyn_bin]);
    }
    for (int dyn_bin = 0; dyn_bin < num_kinematic_bins; ++dyn_bin) {
        if (std::get<0>(resultMh23)[dyn_bin] == 0 ) { continue; }
        // Use average_bin_values instead of bin_center
        aluGraph2.SetPoint(dyn_bin, average_bin_values2[dyn_bin], 
            std::get<0>(resultMh23)[dyn_bin]);
        aluGraph2.SetPointError(dyn_bin, 0, std::get<1>(resultMh23)[dyn_bin]);
    }

    aluGraph1.SetLineColor(kBlue); aluGraph1.SetMarkerColor(kBlue);
    aluGraph1.SetMarkerStyle(20);
    aluGraph1.SetMarkerSize(1.1);

    aluGraph2.SetLineColor(kRed); aluGraph2.SetMarkerColor(kRed);
    aluGraph2.SetMarkerStyle(20);
    aluGraph2.SetMarkerSize(1.1);

    canvas.cd(1);
    aluGraph1.Draw("AP");
    aluGraph1.GetYaxis()->SetRangeUser(-0.15, 0.1);
    aluGraph1.GetYaxis()->SetTitle("F_{LU}^{sin#phi_{#pi^{+}p}} / F_{UU}");
    aluGraph1.GetXaxis()->SetRangeUser(min_val,max_val);
    aluGraph1.GetXaxis()->SetTitle("M_{h(#pi^{+}p)} (GeV)");
    aluGraph1.SetTitle("Q^{2} > 1 (GeV^{2}), W > 2 (GeV), y < 0.75, #it{M}_{X(ep -> e'#pi^{+}#pi^{-}p[X])} < 0.3 (GeV), 0.65 < #it{M}_{h(#pi^{-}X)} (GeV) < 0.91"); 
    aluGraph1.GetXaxis()->SetLabelSize(0.05);  // Increase x-axis label size
    aluGraph1.GetYaxis()->SetLabelSize(0.05);  // Increase y-axis label size
    aluGraph1.GetXaxis()->SetTitleSize(0.07);  // Increase x-axis title size
    aluGraph1.GetYaxis()->SetTitleSize(0.07);  // Increase y-axis title size

    canvas.cd(2);
    aluGraph2.Draw("AP");
    aluGraph2.GetYaxis()->SetRangeUser(-0.15, 0.1);
    aluGraph2.GetYaxis()->SetTitle("F_{LU}^{sin#phi_{#pi^{-}p}} / F_{UU}");
    aluGraph2.GetXaxis()->SetRangeUser(min_val,max_val);
    aluGraph2.GetXaxis()->SetTitle("M_{h(#pi^{-}p)} (GeV)");
    aluGraph2.SetTitle("Q^{2} > 1 (GeV^{2}), W > 2 (GeV), y < 0.75, #it{M}_{X(ep -> e'#pi^{+}#pi^{-}p[X])} < 0.3 (GeV), 0.65 < #it{M}_{h(#pi^{+}X)} < 0.91");  
    aluGraph2.GetXaxis()->SetLabelSize(0.05);  // Increase x-axis label size
    aluGraph2.GetYaxis()->SetLabelSize(0.05);  // Increase y-axis label size
    aluGraph2.GetXaxis()->SetTitleSize(0.07);  // Increase x-axis title size
    aluGraph2.GetYaxis()->SetTitleSize(0.07);  // Increase y-axis title size
    

    // Save the canvas
    canvas.SaveAs("output/BSA_rhominusDeltaplusplus.png");
}

void BSA_rhominusDeltaplusplus(std::string root_file_path) {
    // Start the timer
    auto start_time = std::chrono::high_resolution_clock::now();
    gStyle->SetCanvasColor(0);

    TFile* file = new TFile(root_file_path.c_str(), "READ");

    if (!file->IsOpen()) {
        cout << "Error opening ROOT file (is the location correct?). Exiting." << endl;
    }

    TTree* tree = (TTree*)file->Get("PhysicsEvents");

    if (!tree) {
        cout << "Error getting trees from ROOT file." << endl;
    }

    TTreeReader dataReader(tree); // Create a TTreeReader for the data tree

    createBSAPlot(dataReader, "output");

    file->Close(); delete file;

    // Stop the timer
    auto end_time = std::chrono::high_resolution_clock::now();
    // Calculate the elapsed time in seconds and microseconds
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - 
    start_time).count();
    double seconds = duration / 1e6;
    // Convert to hours, minutes, and seconds
    int hours = static_cast<int>(seconds) / 3600;
    int remaining_time = static_cast<int>(seconds) % 3600;
    int mins = remaining_time / 60;
    int remaining_seconds = remaining_time % 60;
    // Print the elapsed time
    cout << "Time elapsed: ";
    cout << hours << " hours, " << mins << " mins, " << remaining_seconds << " seconds." << endl;
    return 0;
}