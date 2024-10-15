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
    tree->SetBranchStatus("num_neutral", 1);

    // Set branch addresses
    int runnum, num_pos, num_neg, num_neutral;
    tree->SetBranchAddress("runnum", &runnum);
    tree->SetBranchAddress("num_pos", &num_pos);
    tree->SetBranchAddress("num_neg", &num_neg);
    tree->SetBranchAddress("num_neutral", &num_neutral);

    // Define data structures to accumulate sums per run number
    struct Accumulator {
        double sum_pos = 0;
        double sum_pos2 = 0;
        double sum_neg = 0;
        double sum_neg2 = 0;
        double sum_neutral = 0;
        double sum_neutral2 = 0;
        int count = 0;
    };

    std::map<int, Accumulator> run_data;

    // Loop over the tree once
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        Accumulator &acc = run_data[runnum];
        acc.sum_pos += num_pos;
        acc.sum_pos2 += num_pos * num_pos;
        acc.sum_neg += num_neg;
        acc.sum_neg2 += num_neg * num_neg;
        acc.sum_neutral += num_neutral;
        acc.sum_neutral2 += num_neutral * num_neutral;
        acc.count += 1;
    }

    // Prepare vectors for plotting
    std::vector<double> runnum_vals, mean_pos, mean_neg, mean_neutral;
    std::vector<double> err_pos, err_neg, err_neutral;

    // Calculate means and errors for each run number
    for (const auto& kv : run_data) {
        int curr_runnum = kv.first;
        const Accumulator &acc = kv.second;
        int N = acc.count;

        double mean_p = acc.sum_pos / N;
        double mean_n = acc.sum_neg / N;
        double mean_neu = acc.sum_neutral / N;

        // Variance calculation
        double var_p = (acc.sum_pos2 / N) - (mean_p * mean_p);
        double var_n = (acc.sum_neg2 / N) - (mean_n * mean_n);
        double var_neu = (acc.sum_neutral2 / N) - (mean_neu * mean_neu);

        // Statistical uncertainty (standard error)
        double err_p = sqrt(var_p / N);
        double err_n = sqrt(var_n / N);
        double err_neu = sqrt(var_neu / N);

        runnum_vals.push_back(curr_runnum);
        mean_pos.push_back(mean_p);
        mean_neg.push_back(mean_n);
        mean_neutral.push_back(mean_neu);

        err_pos.push_back(err_p);
        err_neg.push_back(err_n);
        err_neutral.push_back(err_neu);
    }

    // Sort data by run number
    std::vector<size_t> indices(runnum_vals.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(),
        [&runnum_vals](size_t i1, size_t i2) { return runnum_vals[i1] < runnum_vals[i2]; });

    // Create sorted vectors
    std::vector<double> runnum_vals_sorted, mean_pos_sorted, mean_neg_sorted, mean_neutral_sorted;
    std::vector<double> err_pos_sorted, err_neg_sorted, err_neutral_sorted;

    for (size_t idx : indices) {
        runnum_vals_sorted.push_back(runnum_vals[idx]);
        mean_pos_sorted.push_back(mean_pos[idx]);
        mean_neg_sorted.push_back(mean_neg[idx]);
        mean_neutral_sorted.push_back(mean_neutral[idx]);
        err_pos_sorted.push_back(err_pos[idx]);
        err_neg_sorted.push_back(err_neg[idx]);
        err_neutral_sorted.push_back(err_neutral[idx]);
    }

    // Create TGraphErrors for plotting
    TGraphErrors *gr_pos = new TGraphErrors(runnum_vals_sorted.size(), &runnum_vals_sorted[0], &mean_pos_sorted[0], 0, &err_pos_sorted[0]);
    TGraphErrors *gr_neg = new TGraphErrors(runnum_vals_sorted.size(), &runnum_vals_sorted[0], &mean_neg_sorted[0], 0, &err_neg_sorted[0]);
    TGraphErrors *gr_neutral = new TGraphErrors(runnum_vals_sorted.size(), &runnum_vals_sorted[0], &mean_neutral_sorted[0], 0, &err_neutral_sorted[0]);

    // Customize graph appearance
    gr_pos->SetMarkerStyle(20);
    gr_pos->SetMarkerColor(kRed);
    gr_pos->SetLineColor(kRed);

    gr_neg->SetMarkerStyle(21);
    gr_neg->SetMarkerColor(kBlue);
    gr_neg->SetLineColor(kBlue);

    gr_neutral->SetMarkerStyle(22);
    gr_neutral->SetMarkerColor(kGreen+2);
    gr_neutral->SetLineColor(kGreen+2);

    // Create and customize canvas
    TCanvas *c1 = new TCanvas("c1", "Multiplicity vs Run Number", 800, 600);
    c1->SetGrid();

    // Draw graphs
    gr_pos->SetTitle("");
    gr_pos->GetXaxis()->SetTitle("runnum");
    gr_pos->GetYaxis()->SetTitle("Multiplicity");
    gr_pos->GetYaxis()->SetRangeUser(0, 2);  // Set y-axis range from 0 to 2
    gr_pos->Draw("AP");
    gr_neg->Draw("P SAME");
    gr_neutral->Draw("P SAME");

    // Add legend
    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gr_pos, "positives", "P");
    leg->AddEntry(gr_neg, "negatives", "P");
    leg->AddEntry(gr_neutral, "neutrals", "P");
    leg->Draw();

    // Save canvas to file
    c1->SaveAs("/home/thayward/nh3_multiplicity.pdf");

    // Clean up
    delete gr_pos;
    delete gr_neg;
    delete gr_neutral;
    delete c1;
    file->Close();
    delete file;
}