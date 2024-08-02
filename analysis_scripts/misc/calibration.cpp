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
#include <sstream>
#include <algorithm>

#include <TGraphErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <sstream>
#include <algorithm>

#include <TGraphErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <sstream>
#include <algorithm>

void plot_htcc_nphe(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    // Arrays to store positive and negative track conditions
    std::vector<int> positive_pids = {-11, 211, 321, 2212};
    std::vector<int> negative_pids = {11, -211, -321, -2212};
    std::vector<int> neutral_pids = {22, 2112};

    // Helper lambda to check if pid is in a vector
    auto is_in = [](int pid, const std::vector<int>& pid_list) {
        return std::find(pid_list.begin(), pid_list.end(), pid) != pid_list.end();
    };

    // Function to create a plot based on track charge
    auto create_plot = [&](const std::string& plot_name, const std::vector<int>& pids) {
        // Restart the TTreeReader to process the data from the beginning
        dataReader.Restart();
        if (mcReader) mcReader->Restart();

        // Set up TTreeReaderValues before calling Next()
        TTreeReaderValue<double> cc_nphe_15(dataReader, "cc_nphe_15");
        TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

        TTreeReaderValue<double>* mc_cc_nphe_15 = nullptr;
        TTreeReaderValue<int>* mc_particle_pid = nullptr;

        if (mcReader) {
            mc_cc_nphe_15 = new TTreeReaderValue<double>(*mcReader, "cc_nphe_15");
            mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
        }

        // Binning for cc_nphe_15
        int nBins = 100;
        double xMin = 0;
        double xMax = 40;

        // Arrays for data
        std::vector<double> dataX(nBins, 0), dataY(nBins, 0);
        std::vector<double> dataEx(nBins, 0), dataEy(nBins, 0);

        // Arrays for MC
        std::vector<double> mcX(nBins, 0), mcY(nBins, 0);
        std::vector<double> mcEx(nBins, 0), mcEy(nBins, 0);

        // Count the number of valid entries
        int dataEntries = 0;
        int mcEntries = 0;

        // Fill the data arrays
        while (dataReader.Next()) {
            double value = *cc_nphe_15;
            int pid = *particle_pid;
            if (value != -9999 && is_in(pid, pids)) {
                int bin = static_cast<int>((value - xMin) / (xMax - xMin) * nBins);
                if (bin >= 0 && bin < nBins) {
                    dataY[bin]++;
                    dataEntries++;
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
                int pid = **mc_particle_pid;
                if (value != -9999 && is_in(pid, pids)) {
                    int bin = static_cast<int>((value - xMin) / (xMax - xMin) * nBins);
                    if (bin >= 0 && bin < nBins) {
                        mcY[bin]++;
                        mcEntries++;
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

        // Determine the maximum Y value for scaling the Y-axis
        double maxY = *std::max_element(dataY.begin(), dataY.end());
        if (mcReader) {
            double maxMCY = *std::max_element(mcY.begin(), mcY.end());
            maxY = std::max(maxY, maxMCY);
        }

        // Create TGraphErrors for data and MC
        TGraphErrors* grData = new 
        	TGraphErrors(nBins, &dataX[0], &dataY[0], &dataEx[0], &dataEy[0]);
        grData->SetMarkerColor(kBlack);
        grData->SetLineColor(kBlack);
        grData->SetTitle(("HTCC nphe - " + plot_name).c_str());
        grData->GetXaxis()->SetTitle("nphe");
        grData->GetYaxis()->SetTitle("normalized counts");

        // Set the Y-axis range to accommodate the maximum value
        grData->SetMaximum(1.15 * maxY);
        grData->GetXaxis()->SetRangeUser(0, 40);

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

        // Create legend with entry counts
        std::ostringstream dataLegend, mcLegend;
        dataLegend << "Data (" << dataEntries << " tracks)";
        TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.9);
        legend->AddEntry(grData, dataLegend.str().c_str(), "pl");
        if (grMC) {
            mcLegend << "MC (" << mcEntries << " tracks)";
            legend->AddEntry(grMC, mcLegend.str().c_str(), "pl");
        }
        legend->Draw();

        // Save the plot
        c.SaveAs(("output/htcc_nphe_" + plot_name + ".png").c_str());

        // Clean up
        delete grData;
        if (grMC) delete grMC;
        delete legend;
        if (mc_cc_nphe_15) delete mc_cc_nphe_15;
        if (mc_particle_pid) delete mc_particle_pid;
    };

    // Create plots for positive and negative tracks
    create_plot("positive", positive_pids);
    create_plot("negative", negative_pids);
}

void plot_ltcc_nphe(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    // Arrays to store positive and negative track conditions
    std::vector<int> positive_pids = {-11, 211, 321, 2212};
    std::vector<int> negative_pids = {11, -211, -321, -2212};
    std::vector<int> neutral_pids = {22, 2112};

    // Helper lambda to check if pid is in a vector
    auto is_in = [](int pid, const std::vector<int>& pid_list) {
        return std::find(pid_list.begin(), pid_list.end(), pid) != pid_list.end();
    };

    // Function to create a plot based on track charge
    auto create_plot = [&](const std::string& plot_name, const std::vector<int>& pids) {
        // Restart the TTreeReader to process the data from the beginning
        dataReader.Restart();
        if (mcReader) mcReader->Restart();

        // Set up TTreeReaderValues before calling Next()
        TTreeReaderValue<double> cc_nphe_16(dataReader, "cc_nphe_16");
        TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

        TTreeReaderValue<double>* mc_cc_nphe_16 = nullptr;
        TTreeReaderValue<int>* mc_particle_pid = nullptr;

        if (mcReader) {
            mc_cc_nphe_16 = new TTreeReaderValue<double>(*mcReader, "cc_nphe_16");
            mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
        }

        // Binning for cc_nphe_16
        int nBins = 100;
        double xMin = 0;
        double xMax = 40;

        // Arrays for data
        std::vector<double> dataX(nBins, 0), dataY(nBins, 0);
        std::vector<double> dataEx(nBins, 0), dataEy(nBins, 0);

        // Arrays for MC
        std::vector<double> mcX(nBins, 0), mcY(nBins, 0);
        std::vector<double> mcEx(nBins, 0), mcEy(nBins, 0);

        // Count the number of valid entries
        int dataEntries = 0;
        int mcEntries = 0;

        // Fill the data arrays
        while (dataReader.Next()) {
            double value = *cc_nphe_16;
            int pid = *particle_pid;
            if (value != -9999 && is_in(pid, pids)) {
                int bin = static_cast<int>((value - xMin) / (xMax - xMin) * nBins);
                if (bin >= 0 && bin < nBins) {
                    dataY[bin]++;
                    dataEntries++;
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
                double value = **mc_cc_nphe_16;
                int pid = **mc_particle_pid;
                if (value != -9999 && is_in(pid, pids)) {
                    int bin = static_cast<int>((value - xMin) / (xMax - xMin) * nBins);
                    if (bin >= 0 && bin < nBins) {
                        mcY[bin]++;
                        mcEntries++;
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

        // Determine the maximum Y value for scaling the Y-axis
        double maxY = *std::max_element(dataY.begin(), dataY.end());
        if (mcReader) {
            double maxMCY = *std::max_element(mcY.begin(), mcY.end());
            maxY = std::max(maxY, maxMCY);
        }

        // Create TGraphErrors for data and MC
        TGraphErrors* grData = new 
        	TGraphErrors(nBins, &dataX[0], &dataY[0], &dataEx[0], &dataEy[0]);
        grData->SetMarkerColor(kBlack);
        grData->SetLineColor(kBlack);
        grData->SetTitle(("LTCC nphe - " + plot_name).c_str());
        grData->GetXaxis()->SetTitle("nphe");
        grData->GetYaxis()->SetTitle("normalized counts");

        // Set the Y-axis range to accommodate the maximum value
        grData->SetMaximum(1.15 * maxY);
        grData->GetXaxis()->SetRangeUser(0, 40);

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

        // Create legend with entry counts
        std::ostringstream dataLegend, mcLegend;
        dataLegend << "Data (" << dataEntries << " tracks)";
        TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.9);
        legend->AddEntry(grData, dataLegend.str().c_str(), "pl");
        if (grMC) {
            mcLegend << "MC (" << mcEntries << " tracks)";
            legend->AddEntry(grMC, mcLegend.str().c_str(), "pl");
        }
        legend->Draw();

        // Save the plot
        c.SaveAs(("output/ltcc_nphe_" + plot_name + ".png").c_str());

        // Clean up
        delete grData;
        if (grMC) delete grMC;
        delete legend;
        if (mc_cc_nphe_16) delete mc_cc_nphe_16;
        if (mc_particle_pid) delete mc_particle_pid;
    };

    // Create plots for positive and negative tracks
    create_plot("positive", positive_pids);
    create_plot("negative", negative_pids);
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

    plot_htcc_nphe(dataReader, mcReader);
    plot_ltcc_nphe(dataReader, mcReader);


    // Close files
    dataFile.Close();
    if (mcFile) {
        mcFile->Close();
        delete mcFile;
    }

    return 0;
}