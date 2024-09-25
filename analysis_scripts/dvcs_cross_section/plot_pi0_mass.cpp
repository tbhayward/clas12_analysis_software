#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <iostream>

// Function to plot pi0 mass from the Mh variable
void plot_pi0_mass(TTreeReader& dataReader1, TTreeReader& dataReader2, TTreeReader& dataReader3,
                   TTreeReader& mcReader1, TTreeReader& mcReader2, TTreeReader& mcReader3,
                   const std::string& outputDir) {
    
    // Set up a 1x3 canvas for the three plots
    TCanvas* canvas = new TCanvas("pi0_mass_canvas", "Invariant Mass of Pi0", 1600, 800);
    canvas->Divide(3, 1);  // 1 row, 3 columns

    // Titles for the three periods
    std::vector<std::string> titles = {"Fa18 Inb", "Fa18 Out", "Sp19 Inb"};

    // Create histograms for each dataset
    std::vector<TH1D*> hist_data;
    std::vector<TH1D*> hist_mc;

    // Bin settings for the invariant mass
    int bins = 50;
    double min_mass = 0.0;
    double max_mass = 0.3;

    // Create histograms for each dataset (data and MC)
    for (int i = 0; i < 3; ++i) {
        hist_data.push_back(new TH1D(("data_" + std::to_string(i)).c_str(), "", bins, min_mass, max_mass));
        hist_mc.push_back(new TH1D(("mc_" + std::to_string(i)).c_str(), "", bins, min_mass, max_mass));
    }

    // TTreeReaders for the Mh (invariant mass) variable
    TTreeReaderValue<double> Mh_data1(dataReader1, "Mh");
    TTreeReaderValue<double> Mh_data2(dataReader2, "Mh");
    TTreeReaderValue<double> Mh_data3(dataReader3, "Mh");

    TTreeReaderValue<double> Mh_mc1(mcReader1, "Mh");
    TTreeReaderValue<double> Mh_mc2(mcReader2, "Mh");
    TTreeReaderValue<double> Mh_mc3(mcReader3, "Mh");

    // Fill histograms for data and MC for each period
    std::vector<TTreeReader*> dataReaders = {&dataReader1, &dataReader2, &dataReader3};
    std::vector<TTreeReader*> mcReaders = {&mcReader1, &mcReader2, &mcReader3};

    for (size_t i = 0; i < 3; ++i) {
        // Fill data histograms
        while (dataReaders[i]->Next()) {
            if (*Mh_data1 > 0.11 && *Mh_data1 < 0.16) {
                hist_data[i]->Fill(*Mh_data1);
                std::cout << *Mh_data1 << std::endl;
            }
        }

        // Fill MC histograms
        while (mcReaders[i]->Next()) {
            if (*Mh_mc1 > 0.11 && *Mh_mc1 < 0.16) {
                hist_mc[i]->Fill(*Mh_mc1);
            }
        }

        // Normalize the histograms
        if (hist_data[i]->Integral() != 0) hist_data[i]->Scale(1.0 / hist_data[i]->Integral());
        if (hist_mc[i]->Integral() != 0) hist_mc[i]->Scale(1.0 / hist_mc[i]->Integral());

        // Set titles and axis labels
        hist_data[i]->SetTitle(titles[i].c_str());
        hist_data[i]->SetXTitle("Invariant Mass (GeV)");
        hist_data[i]->SetYTitle("Normalized Counts");

        hist_mc[i]->SetLineColor(kRed);  // MC in red
        hist_data[i]->SetLineColor(kBlue);  // Data in blue

        // Set the range for the y-axis
        double max_data = hist_data[i]->GetMaximum();
        double max_mc = hist_mc[i]->GetMaximum();
        double y_max = 1.2 * std::max(max_data, max_mc);

        hist_data[i]->GetYaxis()->SetRangeUser(0, y_max);
    }

    // Plot on the canvas
    for (int i = 0; i < 3; ++i) {
        canvas->cd(i + 1);
        hist_data[i]->Draw("HIST");
        hist_mc[i]->Draw("HIST SAME");

        // Add a legend
        TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
        legend->AddEntry(hist_data[i], "Data", "l");
        legend->AddEntry(hist_mc[i], "MC", "l");
        legend->Draw();
    }

    // Save the canvas
    canvas->SaveAs((outputDir + "/pi0_mass/invariant_mass_pi0.png").c_str());

    // Clean up
    delete canvas;
    for (int i = 0; i < 3; ++i) {
        delete hist_data[i];
        delete hist_mc[i];
    }
}