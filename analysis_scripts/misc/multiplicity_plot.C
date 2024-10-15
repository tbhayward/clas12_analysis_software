void multiplicity_plot(const char* filename)
{
    // Open the ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    // Get the TTree named "PhysicsEvents"
    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error: TTree 'PhysicsEvents' not found in file " << filename << std::endl;
        file->Close();
        return;
    }

    // Enable only the necessary branches
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("runnum", 1);
    tree->SetBranchStatus("num_pos", 1);
    tree->SetBranchStatus("num_neg", 1);
    tree->SetBranchStatus("num_neutrals", 1);

    // Set branch addresses
    int runnum, num_pos, num_neg, num_neutrals;
    tree->SetBranchAddress("runnum", &runnum);
    tree->SetBranchAddress("num_pos", &num_pos);
    tree->SetBranchAddress("num_neg", &num_neg);
    tree->SetBranchAddress("num_neutrals", &num_neutrals);

    // Collect unique run numbers
    std::set<int> runnum_set;
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        runnum_set.insert(runnum);
    }

    // Sort run numbers
    std::vector<int> runnums(runnum_set.begin(), runnum_set.end());
    std::sort(runnums.begin(), runnums.end());

    // Vectors to store mean values and errors
    std::vector<double> runnum_vals, mean_pos, mean_neg, mean_neutrals;
    std::vector<double> err_pos, err_neg, err_neutrals;

    // Calculate means and errors for each run number
    for (size_t i = 0; i < runnums.size(); i++) {
        int curr_runnum = runnums[i];
        TString selection = Form("runnum == %d", curr_runnum);

        // Histograms for calculations
        TH1D *h_pos = new TH1D("h_pos", "", 100, -0.5, 99.5);
        TH1D *h_neg = new TH1D("h_neg", "", 100, -0.5, 99.5);
        TH1D *h_neutrals = new TH1D("h_neutrals", "", 100, -0.5, 99.5);

        // Fill histograms
        tree->Draw("num_pos>>h_pos", selection);
        tree->Draw("num_neg>>h_neg", selection);
        tree->Draw("num_neutrals>>h_neutrals", selection);

        // Store mean values and errors
        runnum_vals.push_back(curr_runnum);
        mean_pos.push_back(h_pos->GetMean());
        err_pos.push_back(h_pos->GetMeanError());

        mean_neg.push_back(h_neg->GetMean());
        err_neg.push_back(h_neg->GetMeanError());

        mean_neutrals.push_back(h_neutrals->GetMean());
        err_neutrals.push_back(h_neutrals->GetMeanError());

        // Clean up
        delete h_pos;
        delete h_neg;
        delete h_neutrals;
    }

    // Create TGraphErrors for plotting
    TGraphErrors *gr_pos = new TGraphErrors(runnum_vals.size(), &runnum_vals[0], &mean_pos[0], 0, &err_pos[0]);
    TGraphErrors *gr_neg = new TGraphErrors(runnum_vals.size(), &runnum_vals[0], &mean_neg[0], 0, &err_neg[0]);
    TGraphErrors *gr_neutrals = new TGraphErrors(runnum_vals.size(), &runnum_vals[0], &mean_neutrals[0], 0, &err_neutrals[0]);

    // Customize graph appearance
    gr_pos->SetMarkerStyle(20);
    gr_pos->SetMarkerColor(kRed);
    gr_pos->SetLineColor(kRed);

    gr_neg->SetMarkerStyle(21);
    gr_neg->SetMarkerColor(kBlue);
    gr_neg->SetLineColor(kBlue);

    gr_neutrals->SetMarkerStyle(22);
    gr_neutrals->SetMarkerColor(kGreen+2);
    gr_neutrals->SetLineColor(kGreen+2);

    // Create and customize canvas
    TCanvas *c1 = new TCanvas("c1", "Multiplicity vs Run Number", 800, 600);
    c1->SetGrid();

    // Draw graphs
    gr_pos->SetTitle("");
    gr_pos->GetXaxis()->SetTitle("runnum");
    gr_pos->GetYaxis()->SetTitle("Multiplicity");
    gr_pos->Draw("AP");
    gr_neg->Draw("P SAME");
    gr_neutrals->Draw("P SAME");

    // Add legend
    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gr_pos, "positives", "P");
    leg->AddEntry(gr_neg, "negatives", "P");
    leg->AddEntry(gr_neutrals, "neutrals", "P");
    leg->Draw();

    // Save canvas to file
    c1->SaveAs("/home/thayward/nh3_multiplicity.pdf");

    // Clean up
    delete gr_pos;
    delete gr_neg;
    delete gr_neutrals;
    delete c1;
    file->Close();
    delete file;
}