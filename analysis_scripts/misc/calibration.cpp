#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <iostream>
#include <string>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TF1.h>

void plot_cc_nphe(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    // Set up TTreeReaderValue for cc_nphe_15 in data and MC (if available)
    TTreeReaderValue<double> data_cc_nphe_15(dataReader, "cc_nphe_15");
    TTreeReaderValue<double>* mc_cc_nphe_15 = nullptr;
    if (mcReader) {
        mc_cc_nphe_15 = new TTreeReaderValue<double>(*mcReader, "cc_nphe_15");
    }

    // Binning for cc_nphe_15
    int nBins = 50;
    double xMin = 0;
    double xMax = 50;

    // Arrays for data
    std::vector<double> dataX(nBins, 0), dataY(nBins, 0);
    std::vector<double> dataEx(nBins, 0), dataEy(nBins, 0);

    // Arrays for MC
    std::vector<double> mcX(nBins, 0), mcY(nBins, 0);
    std::vector<double> mcEx(nBins, 0), mcEy(nBins, 0);

    // Fill the data arrays
    while (dataReader.Next()) {
        double value = *data_cc_nphe_15;
        if (value != -9999) {
            int bin = static_cast<int>((value - xMin) / (xMax - xMin) * nBins);
            if (bin >= 0 && bin < nBins) {
                dataY[bin]++;
            }
        }
    }

    // Normalize data counts
    double dataIntegral = 0;
    for (int i = 0; i < nBins; i++) {
        dataIntegral += dataY[i];
    }
    for (int i = 0; i < nBins; i++) {
        dataX[i] = xMin + (i + 0.5) * (xMax - xMin) / nBins;
        dataY[i] /= dataIntegral;
        dataEy[i] = std::sqrt(dataY[i] / dataIntegral);
    }

    // Fill the MC arrays if available
    if (mcReader) {
        while (mcReader->Next()) {
            double value = **mc_cc_nphe_15;
            if (value != -9999) {
                int bin = static_cast<int>((value - xMin) / (xMax - xMin) * nBins);
                if (bin >= 0 && bin < nBins) {
                    mcY[bin]++;
                }
            }
        }

        // Normalize MC counts
        double mcIntegral = 0;
        for (int i = 0; i < nBins; i++) {
            mcIntegral += mcY[i];
        }
        for (int i = 0; i < nBins; i++) {
            mcX[i] = xMin + (i + 0.5) * (xMax - xMin) / nBins;
            mcY[i] /= mcIntegral;
            mcEy[i] = std::sqrt(mcY[i] / mcIntegral);
        }
    }

    // Create TGraphErrors for data and MC
    TGraphErrors* grData = new 
    	TGraphErrors(nBins, &dataX[0], &dataY[0], &dataEx[0], &dataEy[0]);
    grData->SetMarkerColor(kBlack);
    grData->SetLineColor(kBlack);
    grData->SetTitle("cc_nphe_15;nphe;normalized counts");

    TGraphErrors* grMC = nullptr;
    if (mcReader) {
        grMC = new TGraphErrors(nBins, &mcX[0], &mcY[0], &mcEx[0], &mcEy[0]);
        grMC->SetMarkerColor(kRed);
        grMC->SetLineColor(kRed);
    }

    // Draw the plot
    TCanvas c("c", "c", 800, 600);
    grData->Draw("AP");

    if (grMC) {
        grMC->Draw("P SAME");
    }

    // Create legend
    TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(grData, "Data", "pl");
    if (grMC) {
        legend->AddEntry(grMC, "MC", "pl");
    }
    legend->Draw();

    // Save the plot
    c.SaveAs("output/htcc_nphe.png");

    // Clean up
    delete grData;
    if (grMC) delete grMC;
    delete legend;
    if (mc_cc_nphe_15) delete mc_cc_nphe_15;
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << 
        	" <data_file.root> [<mc_file.root>]" << std::endl;
        return 1;
    }

    // Open the first ROOT file
    TFile dataFile(argv[1]);
    if (dataFile.IsZombie()) {
        std::cerr << "Error opening data file!" << std::endl;
        return 1;
    }
    
    // Set up TTreeReader for the first file
    TTreeReader dataReader("PhysicsEvents", &dataFile);
    std::cout << "\nRead in data tree." << std::endl;

    // If a second file is provided, open it and set up a TTreeReader
    TFile* mcFile = nullptr;
    TTreeReader* mcReader = nullptr;
    if (argc == 3) {
        mcFile = new TFile(argv[2]);
        if (mcFile->IsZombie()) {
            std::cerr << "Error opening MC file!" << std::endl;
            return 1;
        }
        mcReader = new TTreeReader("PhysicsEvents", mcFile);
        std::cout << "Read in mc tree." << std::endl;
    }

    //// PLOTS ////
    
    plot_cc_nphe(dataReader, mcReader);


    // Close files
    dataFile.Close();
    if (mcFile) {
        mcFile->Close();
        delete mcFile;
    }

    return 0;
}