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
#include <TH2D.h>

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
        int nBins = 50;
        double xMin = 0;
        double xMax = 20;

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

void plot_ft_xy_energy(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
	gStyle->SetOptStat(0);
    // Restart the TTreeReader to process the data from the beginning
    dataReader.Restart();
    if (mcReader) mcReader->Restart();

    // Set up TTreeReaderValues for ft_x, ft_y, ft_energy, and particle_pid
    TTreeReaderValue<double> ft_x(dataReader, "ft_x");
    TTreeReaderValue<double> ft_y(dataReader, "ft_y");
    TTreeReaderValue<double> ft_energy(dataReader, "ft_energy");
    TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

    TTreeReaderValue<double>* mc_ft_x = nullptr;
    TTreeReaderValue<double>* mc_ft_y = nullptr;
    TTreeReaderValue<double>* mc_ft_energy = nullptr;
    TTreeReaderValue<int>* mc_particle_pid = nullptr;

    if (mcReader) {
        mc_ft_x = new TTreeReaderValue<double>(*mcReader, "ft_x");
        mc_ft_y = new TTreeReaderValue<double>(*mcReader, "ft_y");
        mc_ft_energy = new TTreeReaderValue<double>(*mcReader, "ft_energy");
        mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
    }

    // Define the 2D histogram bins and ranges
    int nBins = 100;
    double xMin = -20;
    double xMax = 20;
    double yMin = -20;
    double yMax = 20;

    // Create histograms for data and MC
    TH2D* h_data_sum = new TH2D("h_data_sum", "Data FT Energy Sum", nBins, xMin, xMax, nBins, yMin, yMax);
    TH2D* h_data_count = new TH2D("h_data_count", "Data FT Count", nBins, xMin, xMax, nBins, yMin, yMax);

    TH2D* h_mc_sum = nullptr;
    TH2D* h_mc_count = nullptr;
    if (mcReader) {
        h_mc_sum = new TH2D("h_mc_sum", "MC FT Energy Sum", nBins, xMin, xMax, nBins, yMin, yMax);
        h_mc_count = new TH2D("h_mc_count", "MC FT Count", nBins, xMin, xMax, nBins, yMin, yMax);
    }

    // Fill the data histograms, applying the cuts
    while (dataReader.Next()) {
        if (*particle_pid == 22 && *ft_x != -9999 && *ft_y != -9999) {
            h_data_sum->Fill(*ft_x, *ft_y, *ft_energy);
            h_data_count->Fill(*ft_x, *ft_y);
        }
    }

    // Fill the MC histograms if available, applying the cuts
    if (mcReader) {
        while (mcReader->Next()) {
            if (**mc_particle_pid == 22 && **mc_ft_x != -9999 && **mc_ft_y != -9999) {
                h_mc_sum->Fill(**mc_ft_x, **mc_ft_y, **mc_ft_energy);
                h_mc_count->Fill(**mc_ft_x, **mc_ft_y);
            }
        }
    }

    // Compute the mean energy for each bin
    TH2D* h_data_mean = new TH2D("h_data_mean", "Data FT Energy Mean", nBins, xMin, xMax, nBins, yMin, yMax);
    h_data_mean->GetXaxis()->SetTitle("x_{FT}");
    h_data_mean->GetYaxis()->SetTitle("y_{FT}");

    if (mcReader) {
        TH2D* h_mc_mean = new TH2D("h_mc_mean", "MC FT Energy Mean", nBins, xMin, xMax, nBins, yMin, yMax);
        h_mc_mean->GetXaxis()->SetTitle("x_{FT}");
        h_mc_mean->GetYaxis()->SetTitle("y_{FT}");

        for (int i = 1; i <= nBins; i++) {
            for (int j = 1; j <= nBins; j++) {
                double count = h_mc_count->GetBinContent(i, j);
                if (count > 0) {
                    h_mc_mean->SetBinContent(i, j, h_mc_sum->GetBinContent(i, j) / count);
                }
            }
        }

        // Compute global mean and standard deviation for MC
        double mc_global_sum = 0, mc_global_count = 0, mc_sum_sq_diff = 0;
        for (int i = 1; i <= nBins; i++) {
            for (int j = 1; j <= nBins; j++) {
                double mean_value = h_mc_mean->GetBinContent(i, j);
                if (mean_value > 0) {
                    mc_global_sum += mean_value;
                    mc_global_count++;
                }
            }
        }
        double mc_global_mean = mc_global_sum / mc_global_count;
        for (int i = 1; i <= nBins; i++) {
            for (int j = 1; j <= nBins; j++) {
                double mean_value = h_mc_mean->GetBinContent(i, j);
                if (mean_value > 0) {
                    mc_sum_sq_diff += TMath::Power(mean_value - mc_global_mean, 2);
                }
            }
        }
        double mc_global_std_dev = TMath::Sqrt(mc_sum_sq_diff / mc_global_count);

        // Draw and save the MC mean energy plot
        TCanvas c_mc("c_mc", "c_mc", 800, 600);
        h_mc_mean->Draw("COLZ");
        TLegend* mc_legend = new TLegend(0.7, 0.8, 0.9, 0.9);
        mc_legend->AddEntry(h_mc_mean, Form("Mean = %.2f GeV", mc_global_mean), "");
        mc_legend->AddEntry(h_mc_mean, Form("Std Dev = %.2f GeV", mc_global_std_dev), "");
        mc_legend->Draw();
        c_mc.SaveAs("output/ft_xy_energy_mc.png");

        // Create masked plot for MC
        TH2D* h_mc_masked = (TH2D*)h_mc_mean->Clone("h_mc_masked");
        TCanvas c_mc_masked("c_mc_masked", "c_mc_masked", 800, 600);
        h_mc_masked->Draw("COLZ");
        for (int i = 1; i <= nBins; i++) {
            for (int j = 1; j <= nBins; j++) {
                double mean_value = h_mc_masked->GetBinContent(i, j);
                if (mean_value < mc_global_mean - 1 * mc_global_std_dev) {
                    TBox* box = new TBox(h_mc_masked->GetXaxis()->GetBinLowEdge(i), h_mc_masked->GetYaxis()->GetBinLowEdge(j),
                                         h_mc_masked->GetXaxis()->GetBinUpEdge(i), h_mc_masked->GetYaxis()->GetBinUpEdge(j));
                    box->SetFillColor(kRed);
                    box->Draw();
                }
            }
        }
        c_mc_masked.SaveAs("output/ft_xy_energy_masked_mc.png");

        // Clean up
        delete mc_legend;
        delete h_mc_masked;
        delete h_mc_sum;
        delete h_mc_count;
        delete h_mc_mean;
    }

    for (int i = 1; i <= nBins; i++) {
        for (int j = 1; j <= nBins; j++) {
            double count = h_data_count->GetBinContent(i, j);
            if (count > 0) {
                h_data_mean->SetBinContent(i, j, h_data_sum->GetBinContent(i, j) / count);
            }
        }
    }

    // Compute global average and standard deviation for Data
    double global_sum = 0, global_count = 0, sum_sq_diff = 0;
    for (int i = 1; i <= nBins; i++) {
        for (int j = 1; j <= nBins; j++) {
            double mean_value = h_data_mean->GetBinContent(i, j);
            if (mean_value > 0) {
                global_sum += mean_value;
                global_count++;
            }
        }
    }
    double global_mean = global_sum / global_count;
    for (int i = 1; i <= nBins; i++) {
        for (int j = 1; j <= nBins; j++) {
            double mean_value = h_data_mean->GetBinContent(i, j);
            if (mean_value > 0) {
                sum_sq_diff += TMath::Power(mean_value - global_mean, 2);
			}
		}
	}
	double global_std_dev = TMath::Sqrt(sum_sq_diff / global_count);
	// Draw and save the data mean energy plot
	TCanvas c_data("c_data", "c_data", 800, 600);
	h_data_mean->Draw("COLZ");
	TLegend* data_legend = new TLegend(0.7, 0.8, 0.9, 0.9);
	data_legend->AddEntry(h_data_mean, Form("Mean = %.2f GeV", global_mean), "");
	data_legend->AddEntry(h_data_mean, Form("Std Dev = %.2f GeV", global_std_dev), "");
	data_legend->Draw();
	c_data.SaveAs("output/ft_xy_energy_data.png");

	// Create masked plot for Data
	TH2D* h_data_masked = (TH2D*)h_data_mean->Clone("h_data_masked");
	TCanvas c_data_masked("c_data_masked", "c_data_masked", 800, 600);
	h_data_masked->Draw("COLZ");
	for (int i = 1; i <= nBins; i++) {
	    for (int j = 1; j <= nBins; j++) {
	        double mean_value = h_data_masked->GetBinContent(i, j);
	        if (mean_value < global_mean - 1 * global_std_dev) {
	            TBox* box = new TBox(h_data_masked->GetXaxis()->GetBinLowEdge(i), h_data_masked->GetYaxis()->GetBinLowEdge(j),
	                                 h_data_masked->GetXaxis()->GetBinUpEdge(i), h_data_masked->GetYaxis()->GetBinUpEdge(j));
	            box->SetFillColor(kRed);
	            box->Draw();
	        }
	    }
	}
	c_data_masked.SaveAs("output/ft_xy_energy_masked_data.png");

	// Clean up
	delete data_legend;
	delete h_data_masked;
	delete h_data_sum;
	delete h_data_count;
	delete h_data_mean;
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

    // plot_htcc_nphe(dataReader, mcReader);
    // plot_ltcc_nphe(dataReader, mcReader);
    // plot_ft_hit_position(dataReader, mcReader);
    plot_ft_xy_energy(dataReader, mcReader);


    // Close files
    dataFile.Close();
    if (mcFile) {
        mcFile->Close();
        delete mcFile;
    }

    return 0;
}