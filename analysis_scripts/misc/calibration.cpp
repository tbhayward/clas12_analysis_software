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
#include <TStyle.h>
#include <TEllipse.h>
#include <TSystem.h>
#include <iostream>
#include <TLine.h> 
#include <TProfile.h>
#include <iostream>
#include <cmath> 
#include <vector>
#include <TLatex.h>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <TPaveText.h>

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
        c.SaveAs(("output/calibration/cc/htcc_nphe_" + plot_name + ".png").c_str());

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
        c.SaveAs(("output/calibration/cc/ltcc_nphe_" + plot_name + ".png").c_str());

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

bool forward_tagger_fiducial(double ft_x, double ft_y) {
    // Compute the radius from the origin
    double radius = sqrt(ft_x * ft_x + ft_y * ft_y);
    
    // Check if the radius is within the fiducial range
    if (radius < 8.5 || radius > 15.5) {
        return false;
    }

    // Define the circle cut parameters: radius and center (x, y)
    std::vector<std::pair<double, std::pair<double, double>>> holes = {
        {1.60, {-8.42,  9.89}},  // circle 1
        {1.60, {-9.89, -5.33}},  // circle 2
        {2.30, {-6.15, -13.00}}, // circle 3
        {2.00, {3.70,  -6.50}}   // circle 4
    };

    // Check if the point falls inside any of the defined holes
    for (const auto& hole : holes) {
        double hole_radius = hole.first;
        double hole_center_x = hole.second.first;
        double hole_center_y = hole.second.second;

        double distance_to_center = sqrt((ft_x - hole_center_x) * (ft_x - hole_center_x) +
                                         (ft_y - hole_center_y) * (ft_y - hole_center_y));

        if (distance_to_center < hole_radius) {
            return false;
        }
    }

    // If all checks are passed, the track is good
    return true;
}

void plot_ft_xy_energy(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    gStyle->SetOptStat(0);

    // Declare and initialize TTreeReaderValue objects before any Next() or Restart() calls
    TTreeReaderValue<double> ft_x(dataReader, "ft_x");
    TTreeReaderValue<double> ft_y(dataReader, "ft_y");
    TTreeReaderValue<double> ft_energy(dataReader, "ft_energy");
    TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

    // Declare pointers for MC TTreeReaderValue objects
    TTreeReaderValue<double>* mc_ft_x = nullptr;
    TTreeReaderValue<double>* mc_ft_y = nullptr;
    TTreeReaderValue<double>* mc_ft_energy = nullptr;
    TTreeReaderValue<int>* mc_particle_pid = nullptr;

    // Initialize MC TTreeReaderValue objects if mcReader is provided
    if (mcReader) {
        mc_ft_x = new TTreeReaderValue<double>(*mcReader, "ft_x");
        mc_ft_y = new TTreeReaderValue<double>(*mcReader, "ft_y");
        mc_ft_energy = new TTreeReaderValue<double>(*mcReader, "ft_energy");
        mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
    }

    // Restart the TTreeReader to process the data from the beginning
    dataReader.Restart();
    if (mcReader) mcReader->Restart();

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

    // Compute the mean energy for each bin (for data)
    TH2D* h_data_mean = new TH2D("h_data_mean", "Data FT Energy Mean", nBins, xMin, xMax, nBins, yMin, yMax);
    h_data_mean->GetXaxis()->SetTitle("x_{FT}");
    h_data_mean->GetYaxis()->SetTitle("y_{FT}");

    for (int i = 1; i <= nBins; i++) {
        for (int j = 1; j <= nBins; j++) {
            double count = h_data_count->GetBinContent(i, j);
            if (count > 0) {
                h_data_mean->SetBinContent(i, j, h_data_sum->GetBinContent(i, j) / count);
            }
        }
    }

    // Compute the mean energy for each bin (for MC)
    TH2D* h_mc_mean = nullptr;
    if (mcReader) {
        h_mc_mean = new TH2D("h_mc_mean", "MC FT Energy Mean", nBins, xMin, xMax, nBins, yMin, yMax);
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
    }

    // Compute global statistics for Data
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

    // Compute global statistics for MC
    double mc_global_sum = 0, mc_global_count = 0, mc_sum_sq_diff = 0;
    double mc_global_mean = 0;
    double mc_global_std_dev = 0;

    if (mcReader) {
        for (int i = 1; i <= nBins; i++) {
            for (int j = 1; j <= nBins; j++) {
                double mean_value = h_mc_mean->GetBinContent(i, j);
                if (mean_value > 0) {
                    mc_global_sum += mean_value;
                    mc_global_count++;
                }
            }
        }
        mc_global_mean = mc_global_sum / mc_global_count;
        for (int i = 1; i <= nBins; i++) {
            for (int j = 1; j <= nBins; j++) {
                double mean_value = h_mc_mean->GetBinContent(i, j);
                if (mean_value > 0) {
                    mc_sum_sq_diff += TMath::Power(mean_value - mc_global_mean, 2);
                }
            }
        }
        mc_global_std_dev = TMath::Sqrt(mc_sum_sq_diff / mc_global_count);
    }

    // Draw and save the data mean energy plot
    TCanvas c_data("c_data", "c_data", 800, 600);
    h_data_mean->Draw("COLZ");
    TLegend* data_legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    data_legend->AddEntry(h_data_mean, Form("Mean = %.2f GeV", global_mean), "");
    data_legend->AddEntry(h_data_mean, Form("Std Dev = %.2f GeV", global_std_dev), "");
    data_legend->Draw();
    c_data.SaveAs("output/calibration/ft/data_ft_xy_energy.png");

    // Draw and save the MC mean energy plot
    TLegend* mc_legend = nullptr;
    if (mcReader) {
        TCanvas c_mc("c_mc", "c_mc", 800, 600);
        h_mc_mean->Draw("COLZ");
        mc_legend = new TLegend(0.7, 0.8, 0.9, 0.9);
        mc_legend->AddEntry(h_mc_mean, Form("Mean = %.2f GeV", mc_global_mean), "");
        mc_legend->AddEntry(h_mc_mean, Form("Std Dev = %.2f GeV", mc_global_std_dev), "");
        mc_legend->Draw();
        c_mc.SaveAs("output/calibration/ft/mc_ft_xy_energy.png");
    }

    // Create and save masked plot for Data
    TH2D* h_data_masked = (TH2D*)h_data_mean->Clone("h_data_masked");
    TCanvas c_data_masked("c_data_masked", "c_data_masked", 800, 600);
    h_data_masked->Draw("COLZ");
    for (int i = 1; i <= nBins; i++) {
        for (int j = 1; j <= nBins; j++) {
            double mean_value = h_data_masked->GetBinContent(i, j);
            if (mean_value < global_mean - 1 * global_std_dev && h_data_mean->GetBinContent(i, j) > 0) {
				TBox* box = new TBox(h_data_masked->GetXaxis()->GetBinLowEdge(i), h_data_masked->GetYaxis()->GetBinLowEdge(j),
				h_data_masked->GetXaxis()->GetBinUpEdge(i), h_data_masked->GetYaxis()->GetBinUpEdge(j));
				box->SetFillColor(kRed);
				box->Draw();
			}
		}
	}

    // Draw the circles representing the holes
	std::vector<std::pair<double, std::pair<double, double>>> holes = {
	    {1.60, {-8.42,  9.89}},  // circle 1
	    {1.60, {-9.89, -5.33}},  // circle 2
	    {2.30, {-6.15, -13.00}},  // circle 3
	    {2.00, {3.70,  -6.50}},   // circle 4
	    {8.5, {0,  0}},   // big circle 1
	    {15.5, {0,  0}},   // big circle 2
	};

	for (size_t idx = 0; idx < holes.size(); ++idx) {
	    double hole_radius = holes[idx].first;
	    double hole_center_x = holes[idx].second.first;
	    double hole_center_y = holes[idx].second.second;

	    TEllipse* circle = new TEllipse(hole_center_x, hole_center_y, hole_radius, hole_radius);
	    circle->SetLineColor(kBlack);
	    circle->SetLineWidth(2);  // Set line width to make it thick
	    circle->SetFillStyle(0);  // No fill color, only the outline
	    
	    // Set dashed line style for the bigger circles
	    if (idx >= 4) {  // Adjust the index based on the order of your circles
	        circle->SetLineStyle(2);  // Dashed line style
	    }

	    circle->Draw("same");
	}

    data_legend->Draw();
    c_data_masked.SaveAs("output/calibration/ft/data_ft_xy_energy_masked.png");

    // Create and save masked plot for MC
	if (mcReader) {
	    TH2D* h_mc_masked = (TH2D*)h_mc_mean->Clone("h_mc_masked");
	    TCanvas c_mc_masked("c_mc_masked", "c_mc_masked", 800, 600);
	    h_mc_masked->Draw("COLZ");
	    for (int i = 1; i <= nBins; i++) {
	        for (int j = 1; j <= nBins; j++) {
	            double mean_value = h_mc_masked->GetBinContent(i, j);
	            if (mean_value < mc_global_mean - 1 * mc_global_std_dev && h_mc_mean->GetBinContent(i, j) > 0) {
	                TBox* box = new TBox(h_mc_masked->GetXaxis()->GetBinLowEdge(i), h_mc_masked->GetYaxis()->GetBinLowEdge(j),
	                                     h_mc_masked->GetXaxis()->GetBinUpEdge(i), h_mc_masked->GetYaxis()->GetBinUpEdge(j));
	                box->SetFillColor(kRed);
	                box->Draw();
	            }
	        }
	    }

	    // Draw the circles representing the holes on MC plot
	    for (size_t idx = 0; idx < holes.size(); ++idx) {
	        double hole_radius = holes[idx].first;
	        double hole_center_x = holes[idx].second.first;
	        double hole_center_y = holes[idx].second.second;

	        TEllipse* circle = new TEllipse(hole_center_x, hole_center_y, hole_radius, hole_radius);
	        circle->SetLineColor(kBlack);
	        circle->SetLineWidth(2);  // Set line width to make it thick
	        circle->SetFillStyle(0);  // No fill color, only the outline

	        // Apply dashed line style only to the last two circles (the larger ones)
	        if (idx >= holes.size() - 2) { // Assuming the last two entries in `holes` are the big circles
	            circle->SetLineStyle(2);  // Dashed line style
	        }

	        circle->Draw("same");
	    }

	    mc_legend->Draw();
	    c_mc_masked.SaveAs("output/calibration/ft/mc_ft_xy_energy_masked.png");
	    delete h_mc_masked;
	}

	// Clean up the dynamically allocated memory
    delete h_data_sum;
    delete h_data_count;
    if (mcReader) {
        delete h_mc_sum;
        delete h_mc_count;
        delete mc_ft_x;
        delete mc_ft_y;
        delete mc_ft_energy;
        delete mc_particle_pid;
    }
}

void plot_ft_hit_position(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {

    // Set up TTreeReaderValues for ft_x, ft_y, and particle_pid
    TTreeReaderValue<double> ft_x(dataReader, "ft_x");
    TTreeReaderValue<double> ft_y(dataReader, "ft_y");
    TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

    TTreeReaderValue<double>* mc_ft_x = nullptr;
    TTreeReaderValue<double>* mc_ft_y = nullptr;
    TTreeReaderValue<int>* mc_particle_pid = nullptr;

    if (mcReader) {
        mc_ft_x = new TTreeReaderValue<double>(*mcReader, "ft_x");
        mc_ft_y = new TTreeReaderValue<double>(*mcReader, "ft_y");
        mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
    }

    // Restart the TTreeReader to process the data from the beginning
    dataReader.Restart();
    if (mcReader) mcReader->Restart();

    // Define the 2D histogram bins and ranges
    int nBins = 100;
    double xMin = -20;
    double xMax = 20;
    double yMin = -20;
    double yMax = 20;

    // Create histograms for data and MC
    TH2D* h_data = new TH2D("h_data", "data FT hit position", nBins, xMin, xMax, nBins, yMin, yMax);
    h_data->GetXaxis()->SetTitle("x_{FT}");
    h_data->GetYaxis()->SetTitle("y_{FT}");

    TH2D* h_mc = nullptr;
    if (mcReader) {
        h_mc = new TH2D("h_mc", "mc FT hit position", nBins, xMin, xMax, nBins, yMin, yMax);
        h_mc->GetXaxis()->SetTitle("x_{FT}");
        h_mc->GetYaxis()->SetTitle("y_{FT}");
    }

    // Create histograms for data and MC with fiducial cuts applied
    TH2D* h_data_cut = new TH2D("h_data_cut", "data FT hit position (cut)", nBins, xMin, xMax, nBins, yMin, yMax);
    h_data_cut->GetXaxis()->SetTitle("x_{FT}");
    h_data_cut->GetYaxis()->SetTitle("y_{FT}");

    TH2D* h_mc_cut = nullptr;
    if (mcReader) {
        h_mc_cut = new TH2D("h_mc_cut", "mc FT hit position (cut)", nBins, xMin, xMax, nBins, yMin, yMax);
        h_mc_cut->GetXaxis()->SetTitle("x_{FT}");
        h_mc_cut->GetYaxis()->SetTitle("y_{FT}");
    }

    // Fill the data histograms, applying the cuts
    while (dataReader.Next()) {
        if (*particle_pid == 22 && *ft_x != -9999 && *ft_y != -9999) {
            h_data->Fill(*ft_x, *ft_y);
            if (forward_tagger_fiducial(*ft_x, *ft_y)) {
                h_data_cut->Fill(*ft_x, *ft_y);
            }
        }
    }

    // Fill the MC histograms if available, applying the cuts
    if (mcReader) {
        while (mcReader->Next()) {
            if (**mc_particle_pid == 22 && **mc_ft_x != -9999 && **mc_ft_y != -9999) {
                h_mc->Fill(**mc_ft_x, **mc_ft_y);
                if (forward_tagger_fiducial(**mc_ft_x, **mc_ft_y)) {
                    h_mc_cut->Fill(**mc_ft_x, **mc_ft_y);
                }
            }
        }
    }

    // Draw and save the original data plot
    TCanvas c_data("c_data", "c_data", 800, 600);
    h_data->Draw("COLZ");
    c_data.SaveAs("output/calibration/ft/data_ft_hit_position.png");

    // Draw and save the original MC plot if available
    if (h_mc) {
        TCanvas c_mc("c_mc", "c_mc", 800, 600);
        h_mc->Draw("COLZ");
        c_mc.SaveAs("output/calibration/ft/mc_ft_hit_position.png");
    }

    // Draw and save the cut data plot
    TCanvas c_data_cut("c_data_cut", "c_data_cut", 800, 600);
    h_data_cut->Draw("COLZ");
    c_data_cut.SaveAs("output/calibration/ft/data_ft_hit_position_cut.png");

    // Draw and save the cut MC plot if available
    if (h_mc_cut) {
        TCanvas c_mc_cut("c_mc_cut", "c_mc_cut", 800, 600);
        h_mc_cut->Draw("COLZ");
        c_mc_cut.SaveAs("output/calibration/ft/mc_ft_hit_position_cut.png");
    }

    // Clean up
    delete h_data;
    delete h_data_cut;
    if (h_mc) delete h_mc;
    if (h_mc_cut) delete h_mc_cut;
    if (mc_ft_x) delete mc_ft_x;
    if (mc_ft_y) delete mc_ft_y;
    if (mc_particle_pid) delete mc_particle_pid;
}

bool pcal_fiducial(double lv_1, double lw_1, double lu_1,
	double lv_4, double lw_4, double lu_4,
	double lv_7, double lw_7, double lu_7,
	int sector, int strictness) {
    // Apply strictness levels for additional cuts
    switch (strictness) {
        case 1:
            if (lw_1 < 9 || lv_1 < 9 || lu_1 < 14) {
                return false;
            }
            break;
        case 2:
            if (lw_1 < 14 || lv_1 < 14 || lu_1 < 29) {
                return false;
            }
            break;
        case 3:
            if ((lw_1 < 19 || lv_1 < 19) || (lu_1 < 39 || lu_1 > 400)) {
                return false;
            }
            break;
        default:
            return false;
    }

    // Specific cuts for each sector in PCal
    // RGA only (not RGC)
    if (sector == 1) {
        if ((lw_1 > 69 && lw_1 < 96) || (lw_1 > 207 && lw_1 < 236)) {
            return false;
        }
    } else if (sector == 2) {
        if ((lv_1 > 95 && lv_1 < 119) || (lu_1 > 108 && lu_1 < 126)) {
            return false;
        }
    } else if (sector == 4) {
        if (lv_1 > 224 && lv_1 < 247) {
            return false;
        }
    } else if (sector == 6) {
        if ((lw_1 > 169 && lw_1 < 198)) {
            return false;
        }
    }

    // Specific cuts for each sector in ECin
    // RGA and RGC
    if (sector == 1) {
    	if (lv_4 > 72 && lv_4 < 94) {
    		return false;
    	}
    }

    // // Specific cuts for each sector in ECout
    // // RGC only so far (not RGA)
    // if (sector == 2) {
    //     if (lw_7 > 68 && lw_7 < 84) {
    //         return false;
    //     }
    // }
    // RGA and RGC
    if (sector == 5) {
    	if (lu_7 > 200 && lu_7 < 220) {
    		return false;
    	}
    }

    // If none of the cuts apply, the track is good
    return true;
}

void plot_pcal_fiducial_determination(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    // Define the 2D histogram bins and ranges
    int nBins_lv_lw_lu = 100;
    int nBins_sf = 40;
    double min = 0;
    double max = 450;
    double lvMin = min;
    double lvMax = max;
    double lwMin = min;
    double lwMax = max;
    double luMin = min;
    double luMax = max;
    double sfMin = 0.15;
    double sfMax = 0.35;

    // Array of particle types (photons and electrons) and their corresponding PIDs
    std::vector<std::tuple<int, std::string>> particle_types = {
        {22, "photon"},
        {11, "electron"}
    };

    // Loop over each particle type
    for (const auto& particle_type : particle_types) {
        int pid = std::get<0>(particle_type);
        std::string particle_name = std::get<1>(particle_type);

        // Restart the TTreeReader to process the data from the beginning
        dataReader.Restart();
        if (mcReader) mcReader->Restart();

        // Declare TTreeReaderValues for data and MC
        TTreeReaderValue<double> cal_energy_1(dataReader, "cal_energy_1");
        TTreeReaderValue<double> cal_energy_4(dataReader, "cal_energy_4");
        TTreeReaderValue<double> cal_energy_7(dataReader, "cal_energy_7");
        TTreeReaderValue<double> p(dataReader, "p");
        TTreeReaderValue<int> cal_sector(dataReader, "cal_sector");
        TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

        TTreeReaderValue<double> cal_lw_1(dataReader, "cal_lw_1");
        TTreeReaderValue<double> cal_lv_1(dataReader, "cal_lv_1");
        TTreeReaderValue<double> cal_lu_1(dataReader, "cal_lu_1");
        TTreeReaderValue<double> cal_lw_4(dataReader, "cal_lw_4");
        TTreeReaderValue<double> cal_lv_4(dataReader, "cal_lv_4");
        TTreeReaderValue<double> cal_lu_4(dataReader, "cal_lu_4");
        TTreeReaderValue<double> cal_lw_7(dataReader, "cal_lw_7");
        TTreeReaderValue<double> cal_lv_7(dataReader, "cal_lv_7");
        TTreeReaderValue<double> cal_lu_7(dataReader, "cal_lu_7");

        TTreeReaderValue<double>* mc_cal_energy_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_energy_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_energy_7 = nullptr;
        TTreeReaderValue<double>* mc_p = nullptr;
        TTreeReaderValue<int>* mc_cal_sector = nullptr;
        TTreeReaderValue<int>* mc_particle_pid = nullptr;

        TTreeReaderValue<double>* mc_cal_lw_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_lv_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_lu_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_lw_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_lv_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_lu_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_lw_7 = nullptr;
        TTreeReaderValue<double>* mc_cal_lv_7 = nullptr;
        TTreeReaderValue<double>* mc_cal_lu_7 = nullptr;

        if (mcReader) {
            mc_cal_energy_1 = new TTreeReaderValue<double>(*mcReader, "cal_energy_1");
            mc_cal_energy_4 = new TTreeReaderValue<double>(*mcReader, "cal_energy_4");
            mc_cal_energy_7 = new TTreeReaderValue<double>(*mcReader, "cal_energy_7");
            mc_p = new TTreeReaderValue<double>(*mcReader, "p");
            mc_cal_sector = new TTreeReaderValue<int>(*mcReader, "cal_sector");
            mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");

            mc_cal_lw_1 = new TTreeReaderValue<double>(*mcReader, "cal_lw_1");
            mc_cal_lv_1 = new TTreeReaderValue<double>(*mcReader, "cal_lv_1");
            mc_cal_lu_1 = new TTreeReaderValue<double>(*mcReader, "cal_lu_1");
            mc_cal_lw_4 = new TTreeReaderValue<double>(*mcReader, "cal_lw_4");
            mc_cal_lv_4 = new TTreeReaderValue<double>(*mcReader, "cal_lv_4");
            mc_cal_lu_4 = new TTreeReaderValue<double>(*mcReader, "cal_lu_4");
            mc_cal_lw_7 = new TTreeReaderValue<double>(*mcReader, "cal_lw_7");
            mc_cal_lv_7 = new TTreeReaderValue<double>(*mcReader, "cal_lv_7");
            mc_cal_lu_7 = new TTreeReaderValue<double>(*mcReader, "cal_lu_7");
        }

        // Create histograms for data and MC for each sector and each combination of lv, lw, lu vs sampling fraction
        TH2D* h_data_lv_lw[6];
        TH2D* h_mc_lv_lw[6] = {nullptr};
        TH2D* h_data_lv_lu[6];
        TH2D* h_mc_lv_lu[6] = {nullptr};
        TH2D* h_data_lw_lu[6];
        TH2D* h_mc_lw_lu[6] = {nullptr};
        TH2D* h_data_sf_lv[6];
        TH2D* h_mc_sf_lv[6] = {nullptr};
        TH2D* h_data_sf_lw[6];
        TH2D* h_mc_sf_lw[6] = {nullptr};
        TH2D* h_data_sf_lu[6];
        TH2D* h_mc_sf_lu[6] = {nullptr};

        for (int sector = 1; sector <= 6; ++sector) {
            std::string title_data_lv_lw = "data PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_lv_lw[sector-1] = new TH2D(("h_data_lv_lw_s" + std::to_string(sector)).c_str(), title_data_lv_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
            h_data_lv_lw[sector-1]->GetXaxis()->SetTitle("lw");
            h_data_lv_lw[sector-1]->GetYaxis()->SetTitle("lv");

            std::string title_data_lv_lu = "data PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_lv_lu[sector-1] = new TH2D(("h_data_lv_lu_s" + std::to_string(sector)).c_str(), title_data_lv_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
            h_data_lv_lu[sector-1]->GetXaxis()->SetTitle("lu");
            h_data_lv_lu[sector-1]->GetYaxis()->SetTitle("lv");

            std::string title_data_lw_lu = "data PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_lw_lu[sector-1] = new TH2D(("h_data_lw_lu_s" + std::to_string(sector)).c_str(), title_data_lw_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
            h_data_lw_lu[sector-1]->GetXaxis()->SetTitle("lu");
            h_data_lw_lu[sector-1]->GetYaxis()->SetTitle("lw");

            std::string title_data_sf_lv = "data PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_sf_lv[sector-1] = new TH2D(("h_data_sf_lv_s" + std::to_string(sector)).c_str(), title_data_sf_lv.c_str(), nBins_lv_lw_lu, lvMin, lvMax, nBins_sf, sfMin, sfMax);
            h_data_sf_lv[sector-1]->GetXaxis()->SetTitle("lv");
            h_data_sf_lv[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

            std::string title_data_sf_lw = "data PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_sf_lw[sector-1] = new TH2D(("h_data_sf_lw_s" + std::to_string(sector)).c_str(), title_data_sf_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_sf, sfMin, sfMax);
            h_data_sf_lw[sector-1]->GetXaxis()->SetTitle("lw");
            h_data_sf_lw[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

            std::string title_data_sf_lu = "data PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_sf_lu[sector-1] = new TH2D(("h_data_sf_lu_s" + std::to_string(sector)).c_str(), title_data_sf_lu.c_str(), nBins_lv_lw_lu, luMin, luMax, nBins_sf, sfMin, sfMax);
            h_data_sf_lu[sector-1]->GetXaxis()->SetTitle("lu");
            h_data_sf_lu[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

            if (mcReader) {
                std::string title_mc_lv_lw = "mc PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_lv_lw[sector-1] = new TH2D(("h_mc_lv_lw_s" + std::to_string(sector)).c_str(), title_mc_lv_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
                h_mc_lv_lw[sector-1]->GetXaxis()->SetTitle("lw");
                h_mc_lv_lw[sector-1]->GetYaxis()->SetTitle("lv");

                std::string title_mc_lv_lu = "mc PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_lv_lu[sector-1] = new TH2D(("h_mc_lv_lu_s" + std::to_string(sector)).c_str(), title_mc_lv_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
                h_mc_lv_lu[sector-1]->GetXaxis()->SetTitle("lu");
                h_mc_lv_lu[sector-1]->GetYaxis()->SetTitle("lv");

                std::string title_mc_lw_lu = "mc PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_lw_lu[sector-1] = new TH2D(("h_mc_lw_lu_s" + std::to_string(sector)).c_str(), title_mc_lw_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
                h_mc_lw_lu[sector-1]->GetXaxis()->SetTitle("lu");
                h_mc_lw_lu[sector-1]->GetYaxis()->SetTitle("lw");

                std::string title_mc_sf_lv = "mc PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_sf_lv[sector-1] = new TH2D(("h_mc_sf_lv_s" + std::to_string(sector)).c_str(), title_mc_sf_lv.c_str(), nBins_lv_lw_lu, lvMin, lvMax, nBins_sf, sfMin, sfMax);
                h_mc_sf_lv[sector-1]->GetXaxis()->SetTitle("lv");
                h_mc_sf_lv[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

                std::string title_mc_sf_lw = "mc PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_sf_lw[sector-1] = new TH2D(("h_mc_sf_lw_s" + std::to_string(sector)).c_str(), title_mc_sf_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_sf, sfMin, sfMax);
                h_mc_sf_lw[sector-1]->GetXaxis()->SetTitle("lw");
                h_mc_sf_lw[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

                std::string title_mc_sf_lu = "mc PCal sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_sf_lu[sector-1] = new TH2D(("h_mc_sf_lu_s" + std::to_string(sector)).c_str(), title_mc_sf_lu.c_str(),nBins_lv_lw_lu, luMin, luMax, nBins_sf, sfMin, sfMax);
                h_mc_sf_lu[sector-1]->GetXaxis()->SetTitle("lu");
                h_mc_sf_lu[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
            }
        }

        // Create TProfile objects for the sampling fraction vs lv, lw, lu
        TProfile* prof_data_sf_lv[6];
        TProfile* prof_data_sf_lw[6];
        TProfile* prof_data_sf_lu[6];

        TProfile* prof_mc_sf_lv[6] = {nullptr};
        TProfile* prof_mc_sf_lw[6] = {nullptr};
        TProfile* prof_mc_sf_lu[6] = {nullptr};

        for (int sector = 0; sector < 6; ++sector) {
            prof_data_sf_lv[sector] = new TProfile(("prof_data_sf_lv_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lv", nBins_lv_lw_lu, lvMin, lvMax, sfMin, sfMax);
            prof_data_sf_lw[sector] = new TProfile(("prof_data_sf_lw_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lw", nBins_lv_lw_lu, lwMin, lwMax, sfMin, sfMax);
            prof_data_sf_lu[sector] = new TProfile(("prof_data_sf_lu_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lu", nBins_lv_lw_lu, luMin, luMax, sfMin, sfMax);

            if (mcReader) {
                prof_mc_sf_lv[sector] = new TProfile(("prof_mc_sf_lv_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lv (MC)", nBins_lv_lw_lu, lvMin, lvMax, sfMin, sfMax);
                prof_mc_sf_lw[sector] = new TProfile(("prof_mc_sf_lw_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lw (MC)", nBins_lv_lw_lu, lwMin, lwMax, sfMin, sfMax);
                prof_mc_sf_lu[sector] = new TProfile(("prof_mc_sf_lu_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lu (MC)", nBins_lv_lw_lu, luMin, luMax, sfMin, sfMax);
            }
        }

        // Fill the histograms and profiles for data
        while (dataReader.Next()) {
            if (*particle_pid == pid && *cal_energy_1 != -9999 && *cal_energy_4 != -9999 && *cal_energy_7 != -9999 && *p != -9999) {
                int sector = *cal_sector;
                if (sector >= 1 && sector <= 6) {
                    double sampling_fraction = (*cal_energy_1 + *cal_energy_4 + *cal_energy_7) / *p;
                    h_data_lv_lw[sector-1]->Fill(*cal_lw_1, *cal_lv_1);
                    h_data_lv_lu[sector-1]->Fill(*cal_lu_1, *cal_lv_1);
                    h_data_lw_lu[sector-1]->Fill(*cal_lu_1, *cal_lw_1);
                    h_data_sf_lv[sector-1]->Fill(*cal_lv_1, sampling_fraction);
                    h_data_sf_lw[sector-1]->Fill(*cal_lw_1, sampling_fraction);
                    h_data_sf_lu[sector-1]->Fill(*cal_lu_1, sampling_fraction);

                    // Fill profiles for comparison plots
                    prof_data_sf_lv[sector-1]->Fill(*cal_lv_1, sampling_fraction);
                    prof_data_sf_lw[sector-1]->Fill(*cal_lw_1, sampling_fraction);
                    prof_data_sf_lu[sector-1]->Fill(*cal_lu_1, sampling_fraction);
                }
            }
        }

        // Fill the histograms and profiles for MC
        if (mcReader) {
            while (mcReader->Next()) {
                if (**mc_particle_pid == pid && **mc_cal_energy_1 != -9999 && **mc_cal_energy_4 != -9999 && **mc_cal_energy_7 != -9999 && **mc_p != -9999) {
                    int sector = **mc_cal_sector;
                    if (sector >= 1 && sector <= 6) {
                        double sampling_fraction = (**mc_cal_energy_1 + **mc_cal_energy_4 + **mc_cal_energy_7) / **mc_p;
                        h_mc_lv_lw[sector-1]->Fill(**mc_cal_lw_1, **mc_cal_lv_1);
                        h_mc_lv_lu[sector-1]->Fill(**mc_cal_lu_1, **mc_cal_lv_1);
                        h_mc_lw_lu[sector-1]->Fill(**mc_cal_lu_1, **mc_cal_lw_1);
                        h_mc_sf_lv[sector-1]->Fill(**mc_cal_lv_1, sampling_fraction);
                        h_mc_sf_lw[sector-1]->Fill(**mc_cal_lw_1, sampling_fraction);
                        h_mc_sf_lu[sector-1]->Fill(**mc_cal_lu_1, sampling_fraction);

                        // Fill profiles for comparison plots
                        prof_mc_sf_lv[sector-1]->Fill(**mc_cal_lv_1, sampling_fraction);
                        prof_mc_sf_lw[sector-1]->Fill(**mc_cal_lw_1, sampling_fraction);
                        prof_mc_sf_lu[sector-1]->Fill(**mc_cal_lu_1, sampling_fraction);
                    }
                }
            }
        }

        // Create TGraphErrors for each sector from TProfile
        TGraphErrors* graph_data_sf_lv[6];
        TGraphErrors* graph_data_sf_lw[6];
        TGraphErrors* graph_data_sf_lu[6];

        TGraphErrors* graph_mc_sf_lv[6] = {nullptr};
        TGraphErrors* graph_mc_sf_lw[6] = {nullptr};
        TGraphErrors* graph_mc_sf_lu[6] = {nullptr};

        for (int sector = 0; sector < 6; ++sector) {
            graph_data_sf_lv[sector] = new TGraphErrors(prof_data_sf_lv[sector]);
            graph_data_sf_lw[sector] = new TGraphErrors(prof_data_sf_lw[sector]);
            graph_data_sf_lu[sector] = new TGraphErrors(prof_data_sf_lu[sector]);

            if (mcReader) {
                graph_mc_sf_lv[sector] = new TGraphErrors(prof_mc_sf_lv[sector]);
                graph_mc_sf_lw[sector] = new TGraphErrors(prof_mc_sf_lw[sector]);
                graph_mc_sf_lu[sector] = new TGraphErrors(prof_mc_sf_lu[sector]);
            }

            // Customize colors
            graph_data_sf_lv[sector]->SetMarkerColor(kBlack);
            graph_data_sf_lv[sector]->SetLineColor(kBlack);
            if (mcReader) {
                graph_mc_sf_lv[sector]->SetMarkerColor(kRed);
                graph_mc_sf_lv[sector]->SetLineColor(kRed);
            }

            graph_data_sf_lw[sector]->SetMarkerColor(kBlack);
            graph_data_sf_lw[sector]->SetLineColor(kBlack);
            if (mcReader) {
                graph_mc_sf_lw[sector]->SetMarkerColor(kRed);
                graph_mc_sf_lw[sector]->SetLineColor(kRed);
            }

            graph_data_sf_lu[sector]->SetMarkerColor(kBlack);
            graph_data_sf_lu[sector]->SetLineColor(kBlack);
            if (mcReader) {
                graph_mc_sf_lu[sector]->SetMarkerColor(kRed);
                graph_mc_sf_lu[sector]->SetLineColor(kRed);
            }
        }

        // Create 2x3 grid plots for comparison of sampling fraction vs lv, lw, lu across all sectors
        TCanvas c_sf_lv_grid("c_sf_lv_grid", ("Sampling Fraction vs lv (all sectors) - " + particle_name).c_str(), 1800, 1200);
        c_sf_lv_grid.Divide(3, 2);
        for (int sector = 1; sector <= 6; ++sector) {
            c_sf_lv_grid.cd(sector);
            gPad->SetLeftMargin(0.20);  // Adjust the left margin to add padding
			graph_data_sf_lv[sector-1]->SetTitle(("Sector " + std::to_string(sector)).c_str());
            graph_data_sf_lv[sector-1]->GetYaxis()->SetRangeUser(0.18, 0.28);
			graph_data_sf_lv[sector-1]->Draw("AP");
			if (mcReader) graph_mc_sf_lv[sector-1]->Draw("P same");
            auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
			legend->AddEntry(graph_data_sf_lv[sector-1], "data", "pl");
			if (mcReader) legend->AddEntry(graph_mc_sf_lv[sector-1], "mc", "pl");
			legend->Draw();
            graph_data_sf_lv[sector-1]->GetXaxis()->SetTitle("lv");
            graph_data_sf_lv[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
        }
        c_sf_lv_grid.SaveAs(("output/calibration/cal/fiducial/pcal/sf_vs_lv_grid_" + particle_name + ".png").c_str());

        TCanvas c_sf_lw_grid("c_sf_lw_grid", ("Sampling Fraction vs lw (all sectors) - " + particle_name).c_str(), 1800, 1200);
        c_sf_lw_grid.Divide(3, 2);
        for (int sector = 1; sector <= 6; ++sector) {
            c_sf_lw_grid.cd(sector);
            gPad->SetLeftMargin(0.20);  // Adjust the left margin to add padding
			graph_data_sf_lw[sector-1]->SetTitle(("Sector " + std::to_string(sector)).c_str());
            graph_data_sf_lw[sector-1]->GetYaxis()->SetRangeUser(0.18, 0.28);
			graph_data_sf_lw[sector-1]->Draw("AP");
			if (mcReader) graph_mc_sf_lw[sector-1]->Draw("P same");
            auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
			legend->AddEntry(graph_data_sf_lw[sector-1], "data", "pl");
			if (mcReader) legend->AddEntry(graph_mc_sf_lv[sector-1], "mc", "pl");
			legend->Draw();
            graph_data_sf_lw[sector-1]->GetXaxis()->SetTitle("lw");
            graph_data_sf_lw[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
        }
        c_sf_lw_grid.SaveAs(("output/calibration/cal/fiducial/pcal/sf_vs_lw_grid_" + particle_name + ".png").c_str());

        TCanvas c_sf_lu_grid("c_sf_lu_grid", ("Sampling Fraction vs lu (all sectors) - " + particle_name).c_str(), 1800, 1200);
        c_sf_lu_grid.Divide(3, 2);
        for (int sector = 1; sector <= 6; ++sector) {
            c_sf_lu_grid.cd(sector);
            gPad->SetLeftMargin(0.20);  // Adjust the left margin to add padding
			graph_data_sf_lu[sector-1]->SetTitle(("Sector " + std::to_string(sector)).c_str());
            graph_data_sf_lu[sector-1]->GetYaxis()->SetRangeUser(0.18, 0.28);
			graph_data_sf_lu[sector-1]->Draw("AP");
			if (mcReader) graph_mc_sf_lu[sector-1]->Draw("P same");
            auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
			legend->AddEntry(graph_data_sf_lu[sector-1], "data", "pl");
			if (mcReader) legend->AddEntry(graph_mc_sf_lu[sector-1], "mc", "pl");
			legend->Draw();
            graph_data_sf_lu[sector-1]->GetXaxis()->SetTitle("lu");
            graph_data_sf_lu[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
        }
        c_sf_lu_grid.SaveAs(("output/calibration/cal/fiducial/pcal/sf_vs_lu_grid_" + particle_name + ".png").c_str());

        // Save the original 2D histograms as before
        TCanvas c_data_lv_lw(("c_data_fiducial_lv_lw_" + particle_name).c_str(), ("Data lv vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_lv_lw.Divide(3, 2);
        c_data_lv_lw.SetLogz(); // Set log scale on the z-axis
        c_data_lv_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        // Save the original 2D histograms as before
        TCanvas c_data_lv_lu(("c_data_fiducial_lv_lu_" + particle_name).c_str(), ("Data lv vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_lv_lu.Divide(3, 2);
        c_data_lv_lu.SetLogz(); // Set log scale on the z-axis
        c_data_lv_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        // Save the original 2D histograms as before
        TCanvas c_data_lw_lu(("c_data_fiducial_lw_lu_" + particle_name).c_str(), ("Data lw vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_lw_lu.Divide(3, 2);
        c_data_lw_lu.SetLogz(); // Set log scale on the z-axis
        c_data_lw_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_lv_lw(("c_mc_fiducial_lv_lw_" + particle_name).c_str(), ("MC lv vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_lv_lw.Divide(3, 2);
        c_mc_lv_lw.SetLogz(); // Set log scale on the z-axis
        c_mc_lv_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_lv_lu(("c_mc_fiducial_lv_lu_" + particle_name).c_str(), ("MC lv vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_lv_lu.Divide(3, 2);
        c_mc_lv_lu.SetLogz(); // Set log scale on the z-axis
        c_mc_lv_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_lw_lu(("c_mc_fiducial_lw_lu_" + particle_name).c_str(), ("MC lw vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_lw_lu.Divide(3, 2);
        c_mc_lw_lu.SetLogz(); // Set log scale on the z-axis
        c_mc_lw_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_data_sf_lv(("c_data_fiducial_sf_lv_" + particle_name).c_str(), ("Data Sampling Fraction vs lv (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_sf_lv.Divide(3, 2);
        c_data_sf_lv.SetLogz(); // Set log scale on the z-axis
        c_data_sf_lv.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_sf_lv(("c_mc_fiducial_sf_lv_" + particle_name).c_str(), ("MC Sampling Fraction vs lv (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_sf_lv.Divide(3, 2);
        c_mc_sf_lv.SetLogz(); // Set log scale on the z-axis
        c_mc_sf_lv.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_data_sf_lw(("c_data_fiducial_sf_lw_" + particle_name).c_str(), ("Data Sampling Fraction vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_sf_lw.Divide(3, 2);
        c_data_sf_lw.SetLogz(); // Set log scale on the z-axis
        c_data_sf_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_sf_lw(("c_mc_fiducial_sf_lw_" + particle_name).c_str(), ("MC Sampling Fraction vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_sf_lw.Divide(3, 2);
        c_mc_sf_lw.SetLogz(); // Set log scale on the z-axis
        c_mc_sf_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_data_sf_lu(("c_data_fiducial_sf_lu_" + particle_name).c_str(), ("Data Sampling Fraction vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_sf_lu.Divide(3, 2);
        c_data_sf_lu.SetLogz(); // Set log scale on the z-axis
        c_data_sf_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_sf_lu(("c_mc_fiducial_sf_lu_" + particle_name).c_str(), ("MC Sampling Fraction vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_sf_lu.Divide(3, 2);
        c_mc_sf_lu.SetLogz(); // Set log scale on the z-axis
        c_mc_sf_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        for (int sector = 1; sector <= 6; ++sector) {
            c_data_lv_lw.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_lv_lw[sector-1]->Draw("COLZ");

            c_data_lv_lu.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_lv_lu[sector-1]->Draw("COLZ");

            c_data_lw_lu.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_lw_lu[sector-1]->Draw("COLZ");

            c_data_sf_lv.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_sf_lv[sector-1]->Draw("COLZ");

            c_data_sf_lw.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_sf_lw[sector-1]->Draw("COLZ");

            c_data_sf_lu.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_sf_lu[sector-1]->Draw("COLZ");

            if (mcReader) {
                c_mc_lv_lw.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_lv_lw[sector-1]->Draw("COLZ");

                c_mc_lv_lu.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_lv_lu[sector-1]->Draw("COLZ");

                c_mc_lw_lu.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_lw_lu[sector-1]->Draw("COLZ");

                c_mc_sf_lv.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_sf_lv[sector-1]->Draw("COLZ");

                c_mc_sf_lw.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_sf_lw[sector-1]->Draw("COLZ");

                c_mc_sf_lu.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_sf_lu[sector-1]->Draw("COLZ");
            }
        }

        // Save the original canvases
        c_data_lv_lw.SaveAs(("output/calibration/cal/fiducial/pcal/data_fiducial_lv_lw_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_lv_lw.SaveAs(("output/calibration/cal/fiducial/pcal/mc_fiducial_lv_lw_" + particle_name + ".png").c_str());

        // Save the original canvases
        c_data_lv_lu.SaveAs(("output/calibration/cal/fiducial/pcal/data_fiducial_lv_lu_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_lv_lu.SaveAs(("output/calibration/cal/fiducial/pcal/mc_fiducial_lv_lu_" + particle_name + ".png").c_str());

        // Save the original canvases
        c_data_lw_lu.SaveAs(("output/calibration/cal/fiducial/pcal/data_fiducial_lw_lu_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_lw_lu.SaveAs(("output/calibration/cal/fiducial/pcal/mc_fiducial_lw_lu_" + particle_name + ".png").c_str());

        c_data_sf_lv.SaveAs(("output/calibration/cal/fiducial/pcal/data_fiducial_sf_lv_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_sf_lv.SaveAs(("output/calibration/cal/fiducial/pcal/mc_fiducial_sf_lv_" + particle_name + ".png").c_str());

        c_data_sf_lw.SaveAs(("output/calibration/cal/fiducial/pcal/data_fiducial_sf_lw_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_sf_lw.SaveAs(("output/calibration/cal/fiducial/pcal/mc_fiducial_sf_lw_" + particle_name + ".png").c_str());

        c_data_sf_lu.SaveAs(("output/calibration/cal/fiducial/pcal/data_fiducial_sf_lu_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_sf_lu.SaveAs(("output/calibration/cal/fiducial/pcal/mc_fiducial_sf_lu_" + particle_name + ".png").c_str());

        // Clean up for this layer and particle type
        for (int sector = 0; sector < 6; ++sector) {
            delete h_data_lv_lw[sector];
            delete h_data_sf_lv[sector];
            delete h_data_sf_lw[sector];
            delete h_data_sf_lu[sector];

            if (mcReader) {
                delete h_mc_lv_lw[sector];
                delete h_mc_sf_lv[sector];


                delete h_mc_sf_lw[sector];
                delete h_mc_sf_lu[sector];
            }
        }

        // Clean up TProfile and TGraphErrors objects
        for (int sector = 0; sector < 6; ++sector) {
            delete prof_data_sf_lv[sector];
            delete prof_data_sf_lw[sector];
            delete prof_data_sf_lu[sector];
            delete graph_data_sf_lv[sector];
            delete graph_data_sf_lw[sector];
            delete graph_data_sf_lu[sector];

            if (mcReader) {
                delete prof_mc_sf_lv[sector];
                delete prof_mc_sf_lw[sector];
                delete prof_mc_sf_lu[sector];
                delete graph_mc_sf_lv[sector];
                delete graph_mc_sf_lw[sector];
                delete graph_mc_sf_lu[sector];
            }
        }

        if (mc_cal_energy_1) delete mc_cal_energy_1;
        if (mc_cal_energy_4) delete mc_cal_energy_4;
        if (mc_cal_energy_7) delete mc_cal_energy_7;
        if (mc_p) delete mc_p;
        if (mc_cal_sector) delete mc_cal_sector;
        if (mc_particle_pid) delete mc_particle_pid;
    }
}

void plot_ecin_fiducial_determination(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    // Define the 2D histogram bins and ranges
    int nBins_lv_lw_lu = 100;
    int nBins_sf = 40;
    double min = 0;
    double max = 450;
    double lvMin = min;
    double lvMax = max;
    double lwMin = min;
    double lwMax = max;
    double luMin = min;
    double luMax = max;
    double sfMin = 0.15;
    double sfMax = 0.35;

    // Array of particle types (photons and electrons) and their corresponding PIDs
    std::vector<std::tuple<int, std::string>> particle_types = {
        {22, "photon"},
        {11, "electron"}
    };

    // Loop over each particle type
    for (const auto& particle_type : particle_types) {
        int pid = std::get<0>(particle_type);
        std::string particle_name = std::get<1>(particle_type);

        // Restart the TTreeReader to process the data from the beginning
        dataReader.Restart();
        if (mcReader) mcReader->Restart();

        // Declare TTreeReaderValues for data and MC
        TTreeReaderValue<double> cal_energy_1(dataReader, "cal_energy_1");
        TTreeReaderValue<double> cal_energy_4(dataReader, "cal_energy_4");
        TTreeReaderValue<double> cal_energy_7(dataReader, "cal_energy_7");
        TTreeReaderValue<double> p(dataReader, "p");
        TTreeReaderValue<int> cal_sector(dataReader, "cal_sector");
        TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

        TTreeReaderValue<double> cal_lw_1(dataReader, "cal_lw_1");
        TTreeReaderValue<double> cal_lv_1(dataReader, "cal_lv_1");
        TTreeReaderValue<double> cal_lu_1(dataReader, "cal_lu_1");
        TTreeReaderValue<double> cal_lw_4(dataReader, "cal_lw_4");
        TTreeReaderValue<double> cal_lv_4(dataReader, "cal_lv_4");
        TTreeReaderValue<double> cal_lu_4(dataReader, "cal_lu_4");
        TTreeReaderValue<double> cal_lw_7(dataReader, "cal_lw_7");
        TTreeReaderValue<double> cal_lv_7(dataReader, "cal_lv_7");
        TTreeReaderValue<double> cal_lu_7(dataReader, "cal_lu_7");

        TTreeReaderValue<double>* mc_cal_energy_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_energy_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_energy_7 = nullptr;
        TTreeReaderValue<double>* mc_p = nullptr;
        TTreeReaderValue<int>* mc_cal_sector = nullptr;
        TTreeReaderValue<int>* mc_particle_pid = nullptr;

        TTreeReaderValue<double>* mc_cal_lw_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_lv_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_lu_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_lw_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_lv_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_lu_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_lw_7 = nullptr;
        TTreeReaderValue<double>* mc_cal_lv_7 = nullptr;
        TTreeReaderValue<double>* mc_cal_lu_7 = nullptr;

        if (mcReader) {
            mc_cal_energy_1 = new TTreeReaderValue<double>(*mcReader, "cal_energy_1");
            mc_cal_energy_4 = new TTreeReaderValue<double>(*mcReader, "cal_energy_4");
            mc_cal_energy_7 = new TTreeReaderValue<double>(*mcReader, "cal_energy_7");
            mc_p = new TTreeReaderValue<double>(*mcReader, "p");
            mc_cal_sector = new TTreeReaderValue<int>(*mcReader, "cal_sector");
            mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");

            mc_cal_lw_1 = new TTreeReaderValue<double>(*mcReader, "cal_lw_1");
            mc_cal_lv_1 = new TTreeReaderValue<double>(*mcReader, "cal_lv_1");
            mc_cal_lu_1 = new TTreeReaderValue<double>(*mcReader, "cal_lu_1");
            mc_cal_lw_4 = new TTreeReaderValue<double>(*mcReader, "cal_lw_4");
            mc_cal_lv_4 = new TTreeReaderValue<double>(*mcReader, "cal_lv_4");
            mc_cal_lu_4 = new TTreeReaderValue<double>(*mcReader, "cal_lu_4");
            mc_cal_lw_7 = new TTreeReaderValue<double>(*mcReader, "cal_lw_7");
            mc_cal_lv_7 = new TTreeReaderValue<double>(*mcReader, "cal_lv_7");
            mc_cal_lu_7 = new TTreeReaderValue<double>(*mcReader, "cal_lu_7");
        }

        // Create histograms for data and MC for each sector and each combination of lv, lw, lu vs sampling fraction
        TH2D* h_data_lv_lw[6];
        TH2D* h_mc_lv_lw[6] = {nullptr};
        TH2D* h_data_lv_lu[6];
        TH2D* h_mc_lv_lu[6] = {nullptr};
        TH2D* h_data_lw_lu[6];
        TH2D* h_mc_lw_lu[6] = {nullptr};
        TH2D* h_data_sf_lv[6];
        TH2D* h_mc_sf_lv[6] = {nullptr};
        TH2D* h_data_sf_lw[6];
        TH2D* h_mc_sf_lw[6] = {nullptr};
        TH2D* h_data_sf_lu[6];
        TH2D* h_mc_sf_lu[6] = {nullptr};

        for (int sector = 1; sector <= 6; ++sector) {
            std::string title_data_lv_lw = "data EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_lv_lw[sector-1] = new TH2D(("h_data_lv_lw_s" + std::to_string(sector)).c_str(), title_data_lv_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
            h_data_lv_lw[sector-1]->GetXaxis()->SetTitle("lw");
            h_data_lv_lw[sector-1]->GetYaxis()->SetTitle("lv");

            std::string title_data_lv_lu = "data EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_lv_lu[sector-1] = new TH2D(("h_data_lv_lu_s" + std::to_string(sector)).c_str(), title_data_lv_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
            h_data_lv_lu[sector-1]->GetXaxis()->SetTitle("lu");
            h_data_lv_lu[sector-1]->GetYaxis()->SetTitle("lv");

            std::string title_data_lw_lu = "data EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_lw_lu[sector-1] = new TH2D(("h_data_lw_lu_s" + std::to_string(sector)).c_str(), title_data_lw_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
            h_data_lw_lu[sector-1]->GetXaxis()->SetTitle("lu");
            h_data_lw_lu[sector-1]->GetYaxis()->SetTitle("lw");

            std::string title_data_sf_lv = "data EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_sf_lv[sector-1] = new TH2D(("h_data_sf_lv_s" + std::to_string(sector)).c_str(), title_data_sf_lv.c_str(), nBins_lv_lw_lu, lvMin, lvMax, nBins_sf, sfMin, sfMax);
            h_data_sf_lv[sector-1]->GetXaxis()->SetTitle("lv");
            h_data_sf_lv[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

            std::string title_data_sf_lw = "data EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_sf_lw[sector-1] = new TH2D(("h_data_sf_lw_s" + std::to_string(sector)).c_str(), title_data_sf_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_sf, sfMin, sfMax);
            h_data_sf_lw[sector-1]->GetXaxis()->SetTitle("lw");
            h_data_sf_lw[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

            std::string title_data_sf_lu = "data EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_sf_lu[sector-1] = new TH2D(("h_data_sf_lu_s" + std::to_string(sector)).c_str(), title_data_sf_lu.c_str(), nBins_lv_lw_lu, luMin, luMax, nBins_sf, sfMin, sfMax);
            h_data_sf_lu[sector-1]->GetXaxis()->SetTitle("lu");
            h_data_sf_lu[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

            if (mcReader) {
                std::string title_mc_lv_lw = "mc EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_lv_lw[sector-1] = new TH2D(("h_mc_lv_lw_s" + std::to_string(sector)).c_str(), title_mc_lv_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
                h_mc_lv_lw[sector-1]->GetXaxis()->SetTitle("lw");
                h_mc_lv_lw[sector-1]->GetYaxis()->SetTitle("lv");

                std::string title_mc_lv_lu = "mc EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_lv_lu[sector-1] = new TH2D(("h_mc_lv_lu_s" + std::to_string(sector)).c_str(), title_mc_lv_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
                h_mc_lv_lu[sector-1]->GetXaxis()->SetTitle("lu");
                h_mc_lv_lu[sector-1]->GetYaxis()->SetTitle("lv");

                std::string title_mc_lw_lu = "mc EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_lw_lu[sector-1] = new TH2D(("h_mc_lw_lu_s" + std::to_string(sector)).c_str(), title_mc_lw_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
                h_mc_lw_lu[sector-1]->GetXaxis()->SetTitle("lu");
                h_mc_lw_lu[sector-1]->GetYaxis()->SetTitle("lw");

                std::string title_mc_sf_lv = "mc EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_sf_lv[sector-1] = new TH2D(("h_mc_sf_lv_s" + std::to_string(sector)).c_str(), title_mc_sf_lv.c_str(), nBins_lv_lw_lu, lvMin, lvMax, nBins_sf, sfMin, sfMax);
                h_mc_sf_lv[sector-1]->GetXaxis()->SetTitle("lv");
                h_mc_sf_lv[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

                std::string title_mc_sf_lw = "mc EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_sf_lw[sector-1] = new TH2D(("h_mc_sf_lw_s" + std::to_string(sector)).c_str(), title_mc_sf_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_sf, sfMin, sfMax);
                h_mc_sf_lw[sector-1]->GetXaxis()->SetTitle("lw");
                h_mc_sf_lw[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

                std::string title_mc_sf_lu = "mc EC_{in} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_sf_lu[sector-1] = new TH2D(("h_mc_sf_lu_s" + std::to_string(sector)).c_str(), title_mc_sf_lu.c_str(),nBins_lv_lw_lu, luMin, luMax, nBins_sf, sfMin, sfMax);
                h_mc_sf_lu[sector-1]->GetXaxis()->SetTitle("lu");
                h_mc_sf_lu[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
            }
        }

        // Create TProfile objects for the sampling fraction vs lv, lw, lu
        TProfile* prof_data_sf_lv[6];
        TProfile* prof_data_sf_lw[6];
        TProfile* prof_data_sf_lu[6];

        TProfile* prof_mc_sf_lv[6] = {nullptr};
        TProfile* prof_mc_sf_lw[6] = {nullptr};
        TProfile* prof_mc_sf_lu[6] = {nullptr};

        for (int sector = 0; sector < 6; ++sector) {
            prof_data_sf_lv[sector] = new TProfile(("prof_data_sf_lv_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lv", nBins_lv_lw_lu, lvMin, lvMax, sfMin, sfMax);
            prof_data_sf_lw[sector] = new TProfile(("prof_data_sf_lw_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lw", nBins_lv_lw_lu, lwMin, lwMax, sfMin, sfMax);
            prof_data_sf_lu[sector] = new TProfile(("prof_data_sf_lu_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lu", nBins_lv_lw_lu, luMin, luMax, sfMin, sfMax);

            if (mcReader) {
                prof_mc_sf_lv[sector] = new TProfile(("prof_mc_sf_lv_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lv (MC)", nBins_lv_lw_lu, lvMin, lvMax, sfMin, sfMax);
                prof_mc_sf_lw[sector] = new TProfile(("prof_mc_sf_lw_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lw (MC)", nBins_lv_lw_lu, lwMin, lwMax, sfMin, sfMax);
                prof_mc_sf_lu[sector] = new TProfile(("prof_mc_sf_lu_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lu (MC)", nBins_lv_lw_lu, luMin, luMax, sfMin, sfMax);
            }
        }

        // Fill the histograms and profiles for data
        while (dataReader.Next()) {
            if (*particle_pid == pid && *cal_energy_1 != -9999 && *cal_energy_4 != -9999 && *cal_energy_7 != -9999 && *p != -9999) {
                int sector = *cal_sector;
                if (sector >= 1 && sector <= 6) {
                    double sampling_fraction = (*cal_energy_1 + *cal_energy_4 + *cal_energy_7) / *p;
                    h_data_lv_lw[sector-1]->Fill(*cal_lw_4, *cal_lv_4);
                    h_data_lv_lu[sector-1]->Fill(*cal_lu_4, *cal_lv_4);
                    h_data_lw_lu[sector-1]->Fill(*cal_lu_4, *cal_lw_4);
                    h_data_sf_lv[sector-1]->Fill(*cal_lv_4, sampling_fraction);
                    h_data_sf_lw[sector-1]->Fill(*cal_lw_4, sampling_fraction);
                    h_data_sf_lu[sector-1]->Fill(*cal_lu_4, sampling_fraction);

                    // Fill profiles for comparison plots
                    prof_data_sf_lv[sector-1]->Fill(*cal_lv_4, sampling_fraction);
                    prof_data_sf_lw[sector-1]->Fill(*cal_lw_4, sampling_fraction);
                    prof_data_sf_lu[sector-1]->Fill(*cal_lu_4, sampling_fraction);
                }
            }
        }

        // Fill the histograms and profiles for MC
        if (mcReader) {
            while (mcReader->Next()) {
                if (**mc_particle_pid == pid && **mc_cal_energy_1 != -9999 && **mc_cal_energy_4 != -9999 && **mc_cal_energy_7 != -9999 && **mc_p != -9999) {
                    int sector = **mc_cal_sector;
                    if (sector >= 1 && sector <= 6) {
                        double sampling_fraction = (**mc_cal_energy_1 + **mc_cal_energy_4 + **mc_cal_energy_7) / **mc_p;
                        h_mc_lv_lw[sector-1]->Fill(**mc_cal_lw_4, **mc_cal_lv_4);
                        h_mc_lv_lu[sector-1]->Fill(**mc_cal_lu_4, **mc_cal_lv_4);
                        h_mc_lw_lu[sector-1]->Fill(**mc_cal_lu_4, **mc_cal_lw_4);
                        h_mc_sf_lv[sector-1]->Fill(**mc_cal_lv_4, sampling_fraction);
                        h_mc_sf_lw[sector-1]->Fill(**mc_cal_lw_4, sampling_fraction);
                        h_mc_sf_lu[sector-1]->Fill(**mc_cal_lu_4, sampling_fraction);

                        // Fill profiles for comparison plots
                        prof_mc_sf_lv[sector-1]->Fill(**mc_cal_lv_4, sampling_fraction);
                        prof_mc_sf_lw[sector-1]->Fill(**mc_cal_lw_4, sampling_fraction);
                        prof_mc_sf_lu[sector-1]->Fill(**mc_cal_lu_4, sampling_fraction);
                    }
                }
            }
        }

        // Create TGraphErrors for each sector from TProfile
        TGraphErrors* graph_data_sf_lv[6];
        TGraphErrors* graph_data_sf_lw[6];
        TGraphErrors* graph_data_sf_lu[6];

        TGraphErrors* graph_mc_sf_lv[6] = {nullptr};
        TGraphErrors* graph_mc_sf_lw[6] = {nullptr};
        TGraphErrors* graph_mc_sf_lu[6] = {nullptr};

        for (int sector = 0; sector < 6; ++sector) {
            graph_data_sf_lv[sector] = new TGraphErrors(prof_data_sf_lv[sector]);
            graph_data_sf_lw[sector] = new TGraphErrors(prof_data_sf_lw[sector]);
            graph_data_sf_lu[sector] = new TGraphErrors(prof_data_sf_lu[sector]);

            if (mcReader) {
                graph_mc_sf_lv[sector] = new TGraphErrors(prof_mc_sf_lv[sector]);
                graph_mc_sf_lw[sector] = new TGraphErrors(prof_mc_sf_lw[sector]);
                graph_mc_sf_lu[sector] = new TGraphErrors(prof_mc_sf_lu[sector]);
            }

            // Customize colors
            graph_data_sf_lv[sector]->SetMarkerColor(kBlack);
            graph_data_sf_lv[sector]->SetLineColor(kBlack);
            if (mcReader) {
                graph_mc_sf_lv[sector]->SetMarkerColor(kRed);
                graph_mc_sf_lv[sector]->SetLineColor(kRed);
            }

            graph_data_sf_lw[sector]->SetMarkerColor(kBlack);
            graph_data_sf_lw[sector]->SetLineColor(kBlack);
            if (mcReader) {
                graph_mc_sf_lw[sector]->SetMarkerColor(kRed);
                graph_mc_sf_lw[sector]->SetLineColor(kRed);
            }

            graph_data_sf_lu[sector]->SetMarkerColor(kBlack);
            graph_data_sf_lu[sector]->SetLineColor(kBlack);
            if (mcReader) {
                graph_mc_sf_lu[sector]->SetMarkerColor(kRed);
                graph_mc_sf_lu[sector]->SetLineColor(kRed);
            }
        }

        // Create 2x3 grid plots for comparison of sampling fraction vs lv, lw, lu across all sectors
        TCanvas c_sf_lv_grid("c_sf_lv_grid", ("Sampling Fraction vs lv (all sectors) - " + particle_name).c_str(), 1800, 1200);
        c_sf_lv_grid.Divide(3, 2);
        for (int sector = 1; sector <= 6; ++sector) {
            c_sf_lv_grid.cd(sector);
            gPad->SetLeftMargin(0.20);  // Adjust the left margin to add padding
			graph_data_sf_lv[sector-1]->SetTitle(("Sector " + std::to_string(sector)).c_str());
            graph_data_sf_lv[sector-1]->GetYaxis()->SetRangeUser(0.18, 0.28);
			graph_data_sf_lv[sector-1]->Draw("AP");
			if (mcReader) graph_mc_sf_lv[sector-1]->Draw("P same");
            auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
			legend->AddEntry(graph_data_sf_lv[sector-1], "data", "pl");
			if (mcReader) legend->AddEntry(graph_mc_sf_lv[sector-1], "mc", "pl");
			legend->Draw();
            graph_data_sf_lv[sector-1]->GetXaxis()->SetTitle("lv");
            graph_data_sf_lv[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
        }
        c_sf_lv_grid.SaveAs(("output/calibration/cal/fiducial/ecin/sf_vs_lv_grid_" + particle_name + ".png").c_str());

        TCanvas c_sf_lw_grid("c_sf_lw_grid", ("Sampling Fraction vs lw (all sectors) - " + particle_name).c_str(), 1800, 1200);
        c_sf_lw_grid.Divide(3, 2);
        for (int sector = 1; sector <= 6; ++sector) {
            c_sf_lw_grid.cd(sector);
            gPad->SetLeftMargin(0.20);  // Adjust the left margin to add padding
			graph_data_sf_lw[sector-1]->SetTitle(("Sector " + std::to_string(sector)).c_str());
            graph_data_sf_lw[sector-1]->GetYaxis()->SetRangeUser(0.18, 0.28);
			graph_data_sf_lw[sector-1]->Draw("AP");
			if (mcReader) graph_mc_sf_lw[sector-1]->Draw("P same");
            auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
			legend->AddEntry(graph_data_sf_lw[sector-1], "data", "pl");
			if (mcReader) legend->AddEntry(graph_mc_sf_lv[sector-1], "mc", "pl");
			legend->Draw();
            graph_data_sf_lw[sector-1]->GetXaxis()->SetTitle("lw");
            graph_data_sf_lw[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
        }
        c_sf_lw_grid.SaveAs(("output/calibration/cal/fiducial/ecin/sf_vs_lw_grid_" + particle_name + ".png").c_str());

        TCanvas c_sf_lu_grid("c_sf_lu_grid", ("Sampling Fraction vs lu (all sectors) - " + particle_name).c_str(), 1800, 1200);
        c_sf_lu_grid.Divide(3, 2);
        for (int sector = 1; sector <= 6; ++sector) {
            c_sf_lu_grid.cd(sector);
            gPad->SetLeftMargin(0.20);  // Adjust the left margin to add padding
			graph_data_sf_lu[sector-1]->SetTitle(("Sector " + std::to_string(sector)).c_str());
            graph_data_sf_lu[sector-1]->GetYaxis()->SetRangeUser(0.18, 0.28);
			graph_data_sf_lu[sector-1]->Draw("AP");
			if (mcReader) graph_mc_sf_lu[sector-1]->Draw("P same");
            auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
			legend->AddEntry(graph_data_sf_lu[sector-1], "data", "pl");
			if (mcReader) legend->AddEntry(graph_mc_sf_lu[sector-1], "mc", "pl");
			legend->Draw();
            graph_data_sf_lu[sector-1]->GetXaxis()->SetTitle("lu");
            graph_data_sf_lu[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
        }
        c_sf_lu_grid.SaveAs(("output/calibration/cal/fiducial/ecin/sf_vs_lu_grid_" + particle_name + ".png").c_str());

        // Save the original 2D histograms as before
        TCanvas c_data_lv_lw(("c_data_fiducial_lv_lw_" + particle_name).c_str(), ("Data lv vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_lv_lw.Divide(3, 2);
        c_data_lv_lw.SetLogz(); // Set log scale on the z-axis
        c_data_lv_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        // Save the original 2D histograms as before
        TCanvas c_data_lv_lu(("c_data_fiducial_lv_lu_" + particle_name).c_str(), ("Data lv vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_lv_lu.Divide(3, 2);
        c_data_lv_lu.SetLogz(); // Set log scale on the z-axis
        c_data_lv_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        // Save the original 2D histograms as before
        TCanvas c_data_lw_lu(("c_data_fiducial_lw_lu_" + particle_name).c_str(), ("Data lw vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_lw_lu.Divide(3, 2);
        c_data_lw_lu.SetLogz(); // Set log scale on the z-axis
        c_data_lw_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_lv_lw(("c_mc_fiducial_lv_lw_" + particle_name).c_str(), ("MC lv vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_lv_lw.Divide(3, 2);
        c_mc_lv_lw.SetLogz(); // Set log scale on the z-axis
        c_mc_lv_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_lv_lu(("c_mc_fiducial_lv_lu_" + particle_name).c_str(), ("MC lv vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_lv_lu.Divide(3, 2);
        c_mc_lv_lu.SetLogz(); // Set log scale on the z-axis
        c_mc_lv_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_lw_lu(("c_mc_fiducial_lw_lu_" + particle_name).c_str(), ("MC lw vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_lw_lu.Divide(3, 2);
        c_mc_lw_lu.SetLogz(); // Set log scale on the z-axis
        c_mc_lw_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_data_sf_lv(("c_data_fiducial_sf_lv_" + particle_name).c_str(), ("Data Sampling Fraction vs lv (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_sf_lv.Divide(3, 2);
        c_data_sf_lv.SetLogz(); // Set log scale on the z-axis
        c_data_sf_lv.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_sf_lv(("c_mc_fiducial_sf_lv_" + particle_name).c_str(), ("MC Sampling Fraction vs lv (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_sf_lv.Divide(3, 2);
        c_mc_sf_lv.SetLogz(); // Set log scale on the z-axis
        c_mc_sf_lv.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_data_sf_lw(("c_data_fiducial_sf_lw_" + particle_name).c_str(), ("Data Sampling Fraction vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_sf_lw.Divide(3, 2);
        c_data_sf_lw.SetLogz(); // Set log scale on the z-axis
        c_data_sf_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_sf_lw(("c_mc_fiducial_sf_lw_" + particle_name).c_str(), ("MC Sampling Fraction vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_sf_lw.Divide(3, 2);
        c_mc_sf_lw.SetLogz(); // Set log scale on the z-axis
        c_mc_sf_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_data_sf_lu(("c_data_fiducial_sf_lu_" + particle_name).c_str(), ("Data Sampling Fraction vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_sf_lu.Divide(3, 2);
        c_data_sf_lu.SetLogz(); // Set log scale on the z-axis
        c_data_sf_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_sf_lu(("c_mc_fiducial_sf_lu_" + particle_name).c_str(), ("MC Sampling Fraction vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_sf_lu.Divide(3, 2);
        c_mc_sf_lu.SetLogz(); // Set log scale on the z-axis
        c_mc_sf_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        for (int sector = 1; sector <= 6; ++sector) {
            c_data_lv_lw.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_lv_lw[sector-1]->Draw("COLZ");

            c_data_lv_lu.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_lv_lu[sector-1]->Draw("COLZ");

            c_data_lw_lu.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_lw_lu[sector-1]->Draw("COLZ");

            c_data_sf_lv.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_sf_lv[sector-1]->Draw("COLZ");

            c_data_sf_lw.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_sf_lw[sector-1]->Draw("COLZ");

            c_data_sf_lu.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_sf_lu[sector-1]->Draw("COLZ");

            if (mcReader) {
                c_mc_lv_lw.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_lv_lw[sector-1]->Draw("COLZ");

                c_mc_lv_lu.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_lv_lu[sector-1]->Draw("COLZ");

                c_mc_lw_lu.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_lw_lu[sector-1]->Draw("COLZ");

                c_mc_sf_lv.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_sf_lv[sector-1]->Draw("COLZ");

                c_mc_sf_lw.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_sf_lw[sector-1]->Draw("COLZ");

                c_mc_sf_lu.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_sf_lu[sector-1]->Draw("COLZ");
            }
        }

        // Save the original canvases
        c_data_lv_lw.SaveAs(("output/calibration/cal/fiducial/ecin/data_fiducial_lv_lw_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_lv_lw.SaveAs(("output/calibration/cal/fiducial/ecin/mc_fiducial_lv_lw_" + particle_name + ".png").c_str());

        // Save the original canvases
        c_data_lv_lu.SaveAs(("output/calibration/cal/fiducial/ecin/data_fiducial_lv_lu_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_lv_lu.SaveAs(("output/calibration/cal/fiducial/ecin/mc_fiducial_lv_lu_" + particle_name + ".png").c_str());

        // Save the original canvases
        c_data_lw_lu.SaveAs(("output/calibration/cal/fiducial/ecin/data_fiducial_lw_lu_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_lw_lu.SaveAs(("output/calibration/cal/fiducial/ecin/mc_fiducial_lw_lu_" + particle_name + ".png").c_str());

        c_data_sf_lv.SaveAs(("output/calibration/cal/fiducial/ecin/data_fiducial_sf_lv_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_sf_lv.SaveAs(("output/calibration/cal/fiducial/ecin/mc_fiducial_sf_lv_" + particle_name + ".png").c_str());

        c_data_sf_lw.SaveAs(("output/calibration/cal/fiducial/ecin/data_fiducial_sf_lw_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_sf_lw.SaveAs(("output/calibration/cal/fiducial/ecin/mc_fiducial_sf_lw_" + particle_name + ".png").c_str());

        c_data_sf_lu.SaveAs(("output/calibration/cal/fiducial/ecin/data_fiducial_sf_lu_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_sf_lu.SaveAs(("output/calibration/cal/fiducial/ecin/mc_fiducial_sf_lu_" + particle_name + ".png").c_str());

        // Clean up for this layer and particle type
        for (int sector = 0; sector < 6; ++sector) {
            delete h_data_lv_lw[sector];
            delete h_data_sf_lv[sector];
            delete h_data_sf_lw[sector];
            delete h_data_sf_lu[sector];

            if (mcReader) {
                delete h_mc_lv_lw[sector];
                delete h_mc_sf_lv[sector];


                delete h_mc_sf_lw[sector];
                delete h_mc_sf_lu[sector];
            }
        }

        // Clean up TProfile and TGraphErrors objects
        for (int sector = 0; sector < 6; ++sector) {
            delete prof_data_sf_lv[sector];
            delete prof_data_sf_lw[sector];
            delete prof_data_sf_lu[sector];
            delete graph_data_sf_lv[sector];
            delete graph_data_sf_lw[sector];
            delete graph_data_sf_lu[sector];

            if (mcReader) {
                delete prof_mc_sf_lv[sector];
                delete prof_mc_sf_lw[sector];
                delete prof_mc_sf_lu[sector];
                delete graph_mc_sf_lv[sector];
                delete graph_mc_sf_lw[sector];
                delete graph_mc_sf_lu[sector];
            }
        }

        if (mc_cal_energy_1) delete mc_cal_energy_1;
        if (mc_cal_energy_4) delete mc_cal_energy_4;
        if (mc_cal_energy_7) delete mc_cal_energy_7;
        if (mc_p) delete mc_p;
        if (mc_cal_sector) delete mc_cal_sector;
        if (mc_particle_pid) delete mc_particle_pid;
    }
}

void plot_ecout_fiducial_determination(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    // Define the 2D histogram bins and ranges
    int nBins_lv_lw_lu = 100;
    int nBins_sf = 40;
    double min = 0;
    double max = 450;
    double lvMin = min;
    double lvMax = max;
    double lwMin = min;
    double lwMax = max;
    double luMin = min;
    double luMax = max;
    double sfMin = 0.15;
    double sfMax = 0.35;

    // Array of particle types (photons and electrons) and their corresponding PIDs
    std::vector<std::tuple<int, std::string>> particle_types = {
        {22, "photon"},
        {11, "electron"}
    };

    // Loop over each particle type
    for (const auto& particle_type : particle_types) {
        int pid = std::get<0>(particle_type);
        std::string particle_name = std::get<1>(particle_type);

        // Restart the TTreeReader to process the data from the beginning
        dataReader.Restart();
        if (mcReader) mcReader->Restart();

        // Declare TTreeReaderValues for data and MC
        TTreeReaderValue<double> cal_energy_1(dataReader, "cal_energy_1");
        TTreeReaderValue<double> cal_energy_4(dataReader, "cal_energy_4");
        TTreeReaderValue<double> cal_energy_7(dataReader, "cal_energy_7");
        TTreeReaderValue<double> p(dataReader, "p");
        TTreeReaderValue<int> cal_sector(dataReader, "cal_sector");
        TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

        TTreeReaderValue<double> cal_lw_1(dataReader, "cal_lw_1");
        TTreeReaderValue<double> cal_lv_1(dataReader, "cal_lv_1");
        TTreeReaderValue<double> cal_lu_1(dataReader, "cal_lu_1");
        TTreeReaderValue<double> cal_lw_4(dataReader, "cal_lw_4");
        TTreeReaderValue<double> cal_lv_4(dataReader, "cal_lv_4");
        TTreeReaderValue<double> cal_lu_4(dataReader, "cal_lu_4");
        TTreeReaderValue<double> cal_lw_7(dataReader, "cal_lw_7");
        TTreeReaderValue<double> cal_lv_7(dataReader, "cal_lv_7");
        TTreeReaderValue<double> cal_lu_7(dataReader, "cal_lu_7");

        TTreeReaderValue<double>* mc_cal_energy_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_energy_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_energy_7 = nullptr;
        TTreeReaderValue<double>* mc_p = nullptr;
        TTreeReaderValue<int>* mc_cal_sector = nullptr;
        TTreeReaderValue<int>* mc_particle_pid = nullptr;

        TTreeReaderValue<double>* mc_cal_lw_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_lv_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_lu_1 = nullptr;
        TTreeReaderValue<double>* mc_cal_lw_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_lv_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_lu_4 = nullptr;
        TTreeReaderValue<double>* mc_cal_lw_7 = nullptr;
        TTreeReaderValue<double>* mc_cal_lv_7 = nullptr;
        TTreeReaderValue<double>* mc_cal_lu_7 = nullptr;

        if (mcReader) {
            mc_cal_energy_1 = new TTreeReaderValue<double>(*mcReader, "cal_energy_1");
            mc_cal_energy_4 = new TTreeReaderValue<double>(*mcReader, "cal_energy_4");
            mc_cal_energy_7 = new TTreeReaderValue<double>(*mcReader, "cal_energy_7");
            mc_p = new TTreeReaderValue<double>(*mcReader, "p");
            mc_cal_sector = new TTreeReaderValue<int>(*mcReader, "cal_sector");
            mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");

            mc_cal_lw_1 = new TTreeReaderValue<double>(*mcReader, "cal_lw_1");
            mc_cal_lv_1 = new TTreeReaderValue<double>(*mcReader, "cal_lv_1");
            mc_cal_lu_1 = new TTreeReaderValue<double>(*mcReader, "cal_lu_1");
            mc_cal_lw_4 = new TTreeReaderValue<double>(*mcReader, "cal_lw_4");
            mc_cal_lv_4 = new TTreeReaderValue<double>(*mcReader, "cal_lv_4");
            mc_cal_lu_4 = new TTreeReaderValue<double>(*mcReader, "cal_lu_4");
            mc_cal_lw_7 = new TTreeReaderValue<double>(*mcReader, "cal_lw_7");
            mc_cal_lv_7 = new TTreeReaderValue<double>(*mcReader, "cal_lv_7");
            mc_cal_lu_7 = new TTreeReaderValue<double>(*mcReader, "cal_lu_7");
        }

        // Create histograms for data and MC for each sector and each combination of lv, lw, lu vs sampling fraction
        TH2D* h_data_lv_lw[6];
        TH2D* h_mc_lv_lw[6] = {nullptr};
        TH2D* h_data_lv_lu[6];
        TH2D* h_mc_lv_lu[6] = {nullptr};
        TH2D* h_data_lw_lu[6];
        TH2D* h_mc_lw_lu[6] = {nullptr};
        TH2D* h_data_sf_lv[6];
        TH2D* h_mc_sf_lv[6] = {nullptr};
        TH2D* h_data_sf_lw[6];
        TH2D* h_mc_sf_lw[6] = {nullptr};
        TH2D* h_data_sf_lu[6];
        TH2D* h_mc_sf_lu[6] = {nullptr};

        for (int sector = 1; sector <= 6; ++sector) {
            std::string title_data_lv_lw = "data EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_lv_lw[sector-1] = new TH2D(("h_data_lv_lw_s" + std::to_string(sector)).c_str(), title_data_lv_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
            h_data_lv_lw[sector-1]->GetXaxis()->SetTitle("lw");
            h_data_lv_lw[sector-1]->GetYaxis()->SetTitle("lv");

            std::string title_data_lv_lu = "data EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_lv_lu[sector-1] = new TH2D(("h_data_lv_lu_s" + std::to_string(sector)).c_str(), title_data_lv_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
            h_data_lv_lu[sector-1]->GetXaxis()->SetTitle("lu");
            h_data_lv_lu[sector-1]->GetYaxis()->SetTitle("lv");

            std::string title_data_lw_lu = "data EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_lw_lu[sector-1] = new TH2D(("h_data_lw_lu_s" + std::to_string(sector)).c_str(), title_data_lw_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
            h_data_lw_lu[sector-1]->GetXaxis()->SetTitle("lu");
            h_data_lw_lu[sector-1]->GetYaxis()->SetTitle("lw");

            std::string title_data_sf_lv = "data EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_sf_lv[sector-1] = new TH2D(("h_data_sf_lv_s" + std::to_string(sector)).c_str(), title_data_sf_lv.c_str(), nBins_lv_lw_lu, lvMin, lvMax, nBins_sf, sfMin, sfMax);
            h_data_sf_lv[sector-1]->GetXaxis()->SetTitle("lv");
            h_data_sf_lv[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

            std::string title_data_sf_lw = "data EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_sf_lw[sector-1] = new TH2D(("h_data_sf_lw_s" + std::to_string(sector)).c_str(), title_data_sf_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_sf, sfMin, sfMax);
            h_data_sf_lw[sector-1]->GetXaxis()->SetTitle("lw");
            h_data_sf_lw[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

            std::string title_data_sf_lu = "data EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
            h_data_sf_lu[sector-1] = new TH2D(("h_data_sf_lu_s" + std::to_string(sector)).c_str(), title_data_sf_lu.c_str(), nBins_lv_lw_lu, luMin, luMax, nBins_sf, sfMin, sfMax);
            h_data_sf_lu[sector-1]->GetXaxis()->SetTitle("lu");
            h_data_sf_lu[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

            if (mcReader) {
                std::string title_mc_lv_lw = "mc EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_lv_lw[sector-1] = new TH2D(("h_mc_lv_lw_s" + std::to_string(sector)).c_str(), title_mc_lv_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
                h_mc_lv_lw[sector-1]->GetXaxis()->SetTitle("lw");
                h_mc_lv_lw[sector-1]->GetYaxis()->SetTitle("lv");

                std::string title_mc_lv_lu = "mc EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_lv_lu[sector-1] = new TH2D(("h_mc_lv_lu_s" + std::to_string(sector)).c_str(), title_mc_lv_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
                h_mc_lv_lu[sector-1]->GetXaxis()->SetTitle("lu");
                h_mc_lv_lu[sector-1]->GetYaxis()->SetTitle("lv");

                std::string title_mc_lw_lu = "mc EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_lw_lu[sector-1] = new TH2D(("h_mc_lw_lu_s" + std::to_string(sector)).c_str(), title_mc_lw_lu.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_lv_lw_lu, lvMin, lvMax);
                h_mc_lw_lu[sector-1]->GetXaxis()->SetTitle("lu");
                h_mc_lw_lu[sector-1]->GetYaxis()->SetTitle("lw");

                std::string title_mc_sf_lv = "mc EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_sf_lv[sector-1] = new TH2D(("h_mc_sf_lv_s" + std::to_string(sector)).c_str(), title_mc_sf_lv.c_str(), nBins_lv_lw_lu, lvMin, lvMax, nBins_sf, sfMin, sfMax);
                h_mc_sf_lv[sector-1]->GetXaxis()->SetTitle("lv");
                h_mc_sf_lv[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

                std::string title_mc_sf_lw = "mc EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_sf_lw[sector-1] = new TH2D(("h_mc_sf_lw_s" + std::to_string(sector)).c_str(), title_mc_sf_lw.c_str(), nBins_lv_lw_lu, lwMin, lwMax, nBins_sf, sfMin, sfMax);
                h_mc_sf_lw[sector-1]->GetXaxis()->SetTitle("lw");
                h_mc_sf_lw[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");

                std::string title_mc_sf_lu = "mc EC_{out} sector " + std::to_string(sector) + " (" + particle_name + ")";
                h_mc_sf_lu[sector-1] = new TH2D(("h_mc_sf_lu_s" + std::to_string(sector)).c_str(), title_mc_sf_lu.c_str(),nBins_lv_lw_lu, luMin, luMax, nBins_sf, sfMin, sfMax);
                h_mc_sf_lu[sector-1]->GetXaxis()->SetTitle("lu");
                h_mc_sf_lu[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
            }
        }

        // Create TProfile objects for the sampling fraction vs lv, lw, lu
        TProfile* prof_data_sf_lv[6];
        TProfile* prof_data_sf_lw[6];
        TProfile* prof_data_sf_lu[6];

        TProfile* prof_mc_sf_lv[6] = {nullptr};
        TProfile* prof_mc_sf_lw[6] = {nullptr};
        TProfile* prof_mc_sf_lu[6] = {nullptr};

        for (int sector = 0; sector < 6; ++sector) {
            prof_data_sf_lv[sector] = new TProfile(("prof_data_sf_lv_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lv", nBins_lv_lw_lu, lvMin, lvMax, sfMin, sfMax);
            prof_data_sf_lw[sector] = new TProfile(("prof_data_sf_lw_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lw", nBins_lv_lw_lu, lwMin, lwMax, sfMin, sfMax);
            prof_data_sf_lu[sector] = new TProfile(("prof_data_sf_lu_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lu", nBins_lv_lw_lu, luMin, luMax, sfMin, sfMax);

            if (mcReader) {
                prof_mc_sf_lv[sector] = new TProfile(("prof_mc_sf_lv_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lv (MC)", nBins_lv_lw_lu, lvMin, lvMax, sfMin, sfMax);
                prof_mc_sf_lw[sector] = new TProfile(("prof_mc_sf_lw_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lw (MC)", nBins_lv_lw_lu, lwMin, lwMax, sfMin, sfMax);
                prof_mc_sf_lu[sector] = new TProfile(("prof_mc_sf_lu_s" + std::to_string(sector+1)).c_str(), "Sampling Fraction vs lu (MC)", nBins_lv_lw_lu, luMin, luMax, sfMin, sfMax);
            }
        }

        // Fill the histograms and profiles for data
        while (dataReader.Next()) {
            if (*particle_pid == pid && *cal_energy_1 != -9999 && *cal_energy_4 != -9999 && *cal_energy_7 != -9999 && *p != -9999) {
                int sector = *cal_sector;
                if (sector >= 1 && sector <= 6) {
                    double sampling_fraction = (*cal_energy_1 + *cal_energy_4 + *cal_energy_7) / *p;
                    h_data_lv_lw[sector-1]->Fill(*cal_lw_7, *cal_lv_7);
                    h_data_lv_lu[sector-1]->Fill(*cal_lu_7, *cal_lv_7);
                    h_data_lw_lu[sector-1]->Fill(*cal_lu_7, *cal_lw_7);
                    h_data_sf_lv[sector-1]->Fill(*cal_lv_7, sampling_fraction);
                    h_data_sf_lw[sector-1]->Fill(*cal_lw_7, sampling_fraction);
                    h_data_sf_lu[sector-1]->Fill(*cal_lu_7, sampling_fraction);

                    // Fill profiles for comparison plots
                    prof_data_sf_lv[sector-1]->Fill(*cal_lv_7, sampling_fraction);
                    prof_data_sf_lw[sector-1]->Fill(*cal_lw_7, sampling_fraction);
                    prof_data_sf_lu[sector-1]->Fill(*cal_lu_7, sampling_fraction);
                }
            }
        }

        // Fill the histograms and profiles for MC
        if (mcReader) {
            while (mcReader->Next()) {
                if (**mc_particle_pid == pid && **mc_cal_energy_1 != -9999 && **mc_cal_energy_4 != -9999 && **mc_cal_energy_7 != -9999 && **mc_p != -9999) {
                    int sector = **mc_cal_sector;
                    if (sector >= 1 && sector <= 6) {
                        double sampling_fraction = (**mc_cal_energy_1 + **mc_cal_energy_4 + **mc_cal_energy_7) / **mc_p;
                        h_mc_lv_lw[sector-1]->Fill(**mc_cal_lw_7, **mc_cal_lv_7);
                        h_mc_lv_lu[sector-1]->Fill(**mc_cal_lu_7, **mc_cal_lv_7);
                        h_mc_lw_lu[sector-1]->Fill(**mc_cal_lu_7, **mc_cal_lw_7);
                        h_mc_sf_lv[sector-1]->Fill(**mc_cal_lv_7, sampling_fraction);
                        h_mc_sf_lw[sector-1]->Fill(**mc_cal_lw_7, sampling_fraction);
                        h_mc_sf_lu[sector-1]->Fill(**mc_cal_lu_7, sampling_fraction);

                        // Fill profiles for comparison plots
                        prof_mc_sf_lv[sector-1]->Fill(**mc_cal_lv_7, sampling_fraction);
                        prof_mc_sf_lw[sector-1]->Fill(**mc_cal_lw_7, sampling_fraction);
                        prof_mc_sf_lu[sector-1]->Fill(**mc_cal_lu_7, sampling_fraction);
                    }
                }
            }
        }

        // Create TGraphErrors for each sector from TProfile
        TGraphErrors* graph_data_sf_lv[6];
        TGraphErrors* graph_data_sf_lw[6];
        TGraphErrors* graph_data_sf_lu[6];

        TGraphErrors* graph_mc_sf_lv[6] = {nullptr};
        TGraphErrors* graph_mc_sf_lw[6] = {nullptr};
        TGraphErrors* graph_mc_sf_lu[6] = {nullptr};

        for (int sector = 0; sector < 6; ++sector) {
            graph_data_sf_lv[sector] = new TGraphErrors(prof_data_sf_lv[sector]);
            graph_data_sf_lw[sector] = new TGraphErrors(prof_data_sf_lw[sector]);
            graph_data_sf_lu[sector] = new TGraphErrors(prof_data_sf_lu[sector]);

            if (mcReader) {
                graph_mc_sf_lv[sector] = new TGraphErrors(prof_mc_sf_lv[sector]);
                graph_mc_sf_lw[sector] = new TGraphErrors(prof_mc_sf_lw[sector]);
                graph_mc_sf_lu[sector] = new TGraphErrors(prof_mc_sf_lu[sector]);
            }

            // Customize colors
            graph_data_sf_lv[sector]->SetMarkerColor(kBlack);
            graph_data_sf_lv[sector]->SetLineColor(kBlack);
            if (mcReader) {
                graph_mc_sf_lv[sector]->SetMarkerColor(kRed);
                graph_mc_sf_lv[sector]->SetLineColor(kRed);
            }

            graph_data_sf_lw[sector]->SetMarkerColor(kBlack);
            graph_data_sf_lw[sector]->SetLineColor(kBlack);
            if (mcReader) {
                graph_mc_sf_lw[sector]->SetMarkerColor(kRed);
                graph_mc_sf_lw[sector]->SetLineColor(kRed);
            }

            graph_data_sf_lu[sector]->SetMarkerColor(kBlack);
            graph_data_sf_lu[sector]->SetLineColor(kBlack);
            if (mcReader) {
                graph_mc_sf_lu[sector]->SetMarkerColor(kRed);
                graph_mc_sf_lu[sector]->SetLineColor(kRed);
            }
        }

        // Create 2x3 grid plots for comparison of sampling fraction vs lv, lw, lu across all sectors
        TCanvas c_sf_lv_grid("c_sf_lv_grid", ("Sampling Fraction vs lv (all sectors) - " + particle_name).c_str(), 1800, 1200);
        c_sf_lv_grid.Divide(3, 2);
        for (int sector = 1; sector <= 6; ++sector) {
            c_sf_lv_grid.cd(sector);
            gPad->SetLeftMargin(0.20);  // Adjust the left margin to add padding
			graph_data_sf_lv[sector-1]->SetTitle(("Sector " + std::to_string(sector)).c_str());
            graph_data_sf_lv[sector-1]->GetYaxis()->SetRangeUser(0.18, 0.28);
			graph_data_sf_lv[sector-1]->Draw("AP");
			if (mcReader) graph_mc_sf_lv[sector-1]->Draw("P same");
            auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
			legend->AddEntry(graph_data_sf_lv[sector-1], "data", "pl");
			if (mcReader) legend->AddEntry(graph_mc_sf_lv[sector-1], "mc", "pl");
			legend->Draw();
            graph_data_sf_lv[sector-1]->GetXaxis()->SetTitle("lv");
            graph_data_sf_lv[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
        }
        c_sf_lv_grid.SaveAs(("output/calibration/cal/fiducial/ecout/sf_vs_lv_grid_" + particle_name + ".png").c_str());

        TCanvas c_sf_lw_grid("c_sf_lw_grid", ("Sampling Fraction vs lw (all sectors) - " + particle_name).c_str(), 1800, 1200);
        c_sf_lw_grid.Divide(3, 2);
        for (int sector = 1; sector <= 6; ++sector) {
            c_sf_lw_grid.cd(sector);
            gPad->SetLeftMargin(0.20);  // Adjust the left margin to add padding
			graph_data_sf_lw[sector-1]->SetTitle(("Sector " + std::to_string(sector)).c_str());
            graph_data_sf_lw[sector-1]->GetYaxis()->SetRangeUser(0.18, 0.28);
			graph_data_sf_lw[sector-1]->Draw("AP");
			if (mcReader) graph_mc_sf_lw[sector-1]->Draw("P same");
            auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
			legend->AddEntry(graph_data_sf_lw[sector-1], "data", "pl");
			if (mcReader) legend->AddEntry(graph_mc_sf_lv[sector-1], "mc", "pl");
			legend->Draw();
            graph_data_sf_lw[sector-1]->GetXaxis()->SetTitle("lw");
            graph_data_sf_lw[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
        }
        c_sf_lw_grid.SaveAs(("output/calibration/cal/fiducial/ecout/sf_vs_lw_grid_" + particle_name + ".png").c_str());

        TCanvas c_sf_lu_grid("c_sf_lu_grid", ("Sampling Fraction vs lu (all sectors) - " + particle_name).c_str(), 1800, 1200);
        c_sf_lu_grid.Divide(3, 2);
        for (int sector = 1; sector <= 6; ++sector) {
            c_sf_lu_grid.cd(sector);
            gPad->SetLeftMargin(0.20);  // Adjust the left margin to add padding
			graph_data_sf_lu[sector-1]->SetTitle(("Sector " + std::to_string(sector)).c_str());
            graph_data_sf_lu[sector-1]->GetYaxis()->SetRangeUser(0.18, 0.28);
			graph_data_sf_lu[sector-1]->Draw("AP");
			if (mcReader) graph_mc_sf_lu[sector-1]->Draw("P same");
            auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
			legend->AddEntry(graph_data_sf_lu[sector-1], "data", "pl");
			if (mcReader) legend->AddEntry(graph_mc_sf_lu[sector-1], "mc", "pl");
			legend->Draw();
            graph_data_sf_lu[sector-1]->GetXaxis()->SetTitle("lu");
            graph_data_sf_lu[sector-1]->GetYaxis()->SetTitle("Sampling Fraction");
        }
        c_sf_lu_grid.SaveAs(("output/calibration/cal/fiducial/ecout/sf_vs_lu_grid_" + particle_name + ".png").c_str());

        // Save the original 2D histograms as before
        TCanvas c_data_lv_lw(("c_data_fiducial_lv_lw_" + particle_name).c_str(), ("Data lv vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_lv_lw.Divide(3, 2);
        c_data_lv_lw.SetLogz(); // Set log scale on the z-axis
        c_data_lv_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        // Save the original 2D histograms as before
        TCanvas c_data_lv_lu(("c_data_fiducial_lv_lu_" + particle_name).c_str(), ("Data lv vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_lv_lu.Divide(3, 2);
        c_data_lv_lu.SetLogz(); // Set log scale on the z-axis
        c_data_lv_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        // Save the original 2D histograms as before
        TCanvas c_data_lw_lu(("c_data_fiducial_lw_lu_" + particle_name).c_str(), ("Data lw vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_lw_lu.Divide(3, 2);
        c_data_lw_lu.SetLogz(); // Set log scale on the z-axis
        c_data_lw_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_lv_lw(("c_mc_fiducial_lv_lw_" + particle_name).c_str(), ("MC lv vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_lv_lw.Divide(3, 2);
        c_mc_lv_lw.SetLogz(); // Set log scale on the z-axis
        c_mc_lv_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_lv_lu(("c_mc_fiducial_lv_lu_" + particle_name).c_str(), ("MC lv vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_lv_lu.Divide(3, 2);
        c_mc_lv_lu.SetLogz(); // Set log scale on the z-axis
        c_mc_lv_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_lw_lu(("c_mc_fiducial_lw_lu_" + particle_name).c_str(), ("MC lw vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_lw_lu.Divide(3, 2);
        c_mc_lw_lu.SetLogz(); // Set log scale on the z-axis
        c_mc_lw_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_data_sf_lv(("c_data_fiducial_sf_lv_" + particle_name).c_str(), ("Data Sampling Fraction vs lv (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_sf_lv.Divide(3, 2);
        c_data_sf_lv.SetLogz(); // Set log scale on the z-axis
        c_data_sf_lv.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_sf_lv(("c_mc_fiducial_sf_lv_" + particle_name).c_str(), ("MC Sampling Fraction vs lv (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_sf_lv.Divide(3, 2);
        c_mc_sf_lv.SetLogz(); // Set log scale on the z-axis
        c_mc_sf_lv.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_data_sf_lw(("c_data_fiducial_sf_lw_" + particle_name).c_str(), ("Data Sampling Fraction vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_sf_lw.Divide(3, 2);
        c_data_sf_lw.SetLogz(); // Set log scale on the z-axis
        c_data_sf_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_sf_lw(("c_mc_fiducial_sf_lw_" + particle_name).c_str(), ("MC Sampling Fraction vs lw (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_sf_lw.Divide(3, 2);
        c_mc_sf_lw.SetLogz(); // Set log scale on the z-axis
        c_mc_sf_lw.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_data_sf_lu(("c_data_fiducial_sf_lu_" + particle_name).c_str(), ("Data Sampling Fraction vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_data_sf_lu.Divide(3, 2);
        c_data_sf_lu.SetLogz(); // Set log scale on the z-axis
        c_data_sf_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        TCanvas c_mc_sf_lu(("c_mc_fiducial_sf_lu_" + particle_name).c_str(), ("MC Sampling Fraction vs lu (" + particle_name + ")").c_str(), 1800, 1200);
        c_mc_sf_lu.Divide(3, 2);
        c_mc_sf_lu.SetLogz(); // Set log scale on the z-axis
        c_mc_sf_lu.SetMargin(0.15, 0.15, 0.15, 0.15); // Add padding

        for (int sector = 1; sector <= 6; ++sector) {
            c_data_lv_lw.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_lv_lw[sector-1]->Draw("COLZ");

            c_data_lv_lu.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_lv_lu[sector-1]->Draw("COLZ");

            c_data_lw_lu.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_lw_lu[sector-1]->Draw("COLZ");

            c_data_sf_lv.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_sf_lv[sector-1]->Draw("COLZ");

            c_data_sf_lw.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_sf_lw[sector-1]->Draw("COLZ");

            c_data_sf_lu.cd(sector);
            gPad->SetLogz(); // Set log scale on the z-axis
            gPad->SetLeftMargin(0.15); // Adjust the left margin
            gPad->SetRightMargin(0.15); // Adjust the right margin
            h_data_sf_lu[sector-1]->Draw("COLZ");

            if (mcReader) {
                c_mc_lv_lw.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_lv_lw[sector-1]->Draw("COLZ");

                c_mc_lv_lu.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_lv_lu[sector-1]->Draw("COLZ");

                c_mc_lw_lu.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_lw_lu[sector-1]->Draw("COLZ");

                c_mc_sf_lv.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_sf_lv[sector-1]->Draw("COLZ");

                c_mc_sf_lw.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_sf_lw[sector-1]->Draw("COLZ");

                c_mc_sf_lu.cd(sector);
                gPad->SetLogz(); // Set log scale on the z-axis
                gPad->SetLeftMargin(0.15); // Adjust the left margin
                gPad->SetRightMargin(0.15); // Adjust the right margin
                h_mc_sf_lu[sector-1]->Draw("COLZ");
            }
        }

        // Save the original canvases
        c_data_lv_lw.SaveAs(("output/calibration/cal/fiducial/ecout/data_fiducial_lv_lw_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_lv_lw.SaveAs(("output/calibration/cal/fiducial/ecout//mc_fiducial_lv_lw_" + particle_name + ".png").c_str());

        // Save the original canvases
        c_data_lv_lu.SaveAs(("output/calibration/cal/fiducial/ecout/data_fiducial_lv_lu_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_lv_lu.SaveAs(("output/calibration/cal/fiducial/ecout//mc_fiducial_lv_lu_" + particle_name + ".png").c_str());

        // Save the original canvases
        c_data_lw_lu.SaveAs(("output/calibration/cal/fiducial/ecout/data_fiducial_lw_lu_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_lw_lu.SaveAs(("output/calibration/cal/fiducial/ecout/mc_fiducial_lw_lu_" + particle_name + ".png").c_str());

        c_data_sf_lv.SaveAs(("output/calibration/cal/fiducial/ecout/data_fiducial_sf_lv_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_sf_lv.SaveAs(("output/calibration/cal/fiducial/ecout/mc_fiducial_sf_lv_" + particle_name + ".png").c_str());

        c_data_sf_lw.SaveAs(("output/calibration/cal/fiducial/ecout/data_fiducial_sf_lw_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_sf_lw.SaveAs(("output/calibration/cal/fiducial/ecout/mc_fiducial_sf_lw_" + particle_name + ".png").c_str());

        c_data_sf_lu.SaveAs(("output/calibration/cal/fiducial/ecout/data_fiducial_sf_lu_" + particle_name + ".png").c_str());
        if (mcReader) c_mc_sf_lu.SaveAs(("output/calibration/cal/fiducial/ecout/mc_fiducial_sf_lu_" + particle_name + ".png").c_str());

        // Clean up for this layer and particle type
        for (int sector = 0; sector < 6; ++sector) {
            delete h_data_lv_lw[sector];
            delete h_data_sf_lv[sector];
            delete h_data_sf_lw[sector];
            delete h_data_sf_lu[sector];

            if (mcReader) {
                delete h_mc_lv_lw[sector];
                delete h_mc_sf_lv[sector];


                delete h_mc_sf_lw[sector];
                delete h_mc_sf_lu[sector];
            }
        }

        // Clean up TProfile and TGraphErrors objects
        for (int sector = 0; sector < 6; ++sector) {
            delete prof_data_sf_lv[sector];
            delete prof_data_sf_lw[sector];
            delete prof_data_sf_lu[sector];
            delete graph_data_sf_lv[sector];
            delete graph_data_sf_lw[sector];
            delete graph_data_sf_lu[sector];

            if (mcReader) {
                delete prof_mc_sf_lv[sector];
                delete prof_mc_sf_lw[sector];
                delete prof_mc_sf_lu[sector];
                delete graph_mc_sf_lv[sector];
                delete graph_mc_sf_lw[sector];
                delete graph_mc_sf_lu[sector];
            }
        }

        if (mc_cal_energy_1) delete mc_cal_energy_1;
        if (mc_cal_energy_4) delete mc_cal_energy_4;
        if (mc_cal_energy_7) delete mc_cal_energy_7;
        if (mc_p) delete mc_p;
        if (mc_cal_sector) delete mc_cal_sector;
        if (mc_particle_pid) delete mc_particle_pid;
    }
}

void plot_cal_hit_position(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    // Define the 2D histogram bins and ranges
    int nBins = 100;
    double xMin = -450;
    double xMax = 450;
    double yMin = -450;
    double yMax = 450;

    // Array of layers and their corresponding names
    std::vector<std::tuple<std::string, std::string, std::string>> layers = {
        {"cal_x_1", "cal_y_1", "PCal"},
        {"cal_x_4", "cal_y_4", "EC_{in}"},
        {"cal_x_7", "cal_y_7", "EC_{out}"}
    };

    // Array of particle types (photons and electrons) and their corresponding PIDs
    std::vector<std::tuple<int, std::string>> particle_types = {
        {22, "photon"},
        {11, "electron"}
    };

    // Declare TTreeReaderValues for the layer 1 cuts (lv, lw, lu)
    TTreeReaderValue<double> cal_lv_1(dataReader, "cal_lv_1");
    TTreeReaderValue<double> cal_lw_1(dataReader, "cal_lw_1");
    TTreeReaderValue<double> cal_lu_1(dataReader, "cal_lu_1");
    TTreeReaderValue<double> cal_lv_4(dataReader, "cal_lv_4");
    TTreeReaderValue<double> cal_lw_4(dataReader, "cal_lw_4");
    TTreeReaderValue<double> cal_lu_4(dataReader, "cal_lu_4");
    TTreeReaderValue<double> cal_lv_7(dataReader, "cal_lv_7");
    TTreeReaderValue<double> cal_lw_7(dataReader, "cal_lw_7");
    TTreeReaderValue<double> cal_lu_7(dataReader, "cal_lu_7");
    TTreeReaderValue<int> cal_sector(dataReader, "cal_sector");

    TTreeReaderValue<double>* mc_cal_lv_1 = nullptr;
    TTreeReaderValue<double>* mc_cal_lw_1 = nullptr;
    TTreeReaderValue<double>* mc_cal_lu_1 = nullptr;
    TTreeReaderValue<double>* mc_cal_lv_4 = nullptr;
    TTreeReaderValue<double>* mc_cal_lw_4 = nullptr;
    TTreeReaderValue<double>* mc_cal_lu_4 = nullptr;
    TTreeReaderValue<double>* mc_cal_lv_7 = nullptr;
    TTreeReaderValue<double>* mc_cal_lw_7 = nullptr;
    TTreeReaderValue<double>* mc_cal_lu_7 = nullptr;
    TTreeReaderValue<int>* mc_cal_sector = nullptr;

    if (mcReader) {
        mc_cal_lv_1 = new TTreeReaderValue<double>(*mcReader, "cal_lv_1");
        mc_cal_lw_1 = new TTreeReaderValue<double>(*mcReader, "cal_lw_1");
        mc_cal_lu_1 = new TTreeReaderValue<double>(*mcReader, "cal_lu_1");
        mc_cal_lv_4 = new TTreeReaderValue<double>(*mcReader, "cal_lv_4");
        mc_cal_lw_4 = new TTreeReaderValue<double>(*mcReader, "cal_lw_4");
        mc_cal_lu_4 = new TTreeReaderValue<double>(*mcReader, "cal_lu_4");
        mc_cal_lv_7 = new TTreeReaderValue<double>(*mcReader, "cal_lv_7");
        mc_cal_lw_7 = new TTreeReaderValue<double>(*mcReader, "cal_lw_7");
        mc_cal_lu_7 = new TTreeReaderValue<double>(*mcReader, "cal_lu_7");
        mc_cal_sector = new TTreeReaderValue<int>(*mcReader, "cal_sector");
    }

    // Loop over each particle type
    for (const auto& particle_type : particle_types) {
        int pid = std::get<0>(particle_type);
        std::string particle_name = std::get<1>(particle_type);

        // Loop over each layer
        for (const auto& layer : layers) {
            std::string x_branch = std::get<0>(layer);
            std::string y_branch = std::get<1>(layer);
            std::string layer_name = std::get<2>(layer);

            // Restart the TTreeReader to process the data from the beginning
            dataReader.Restart();
            if (mcReader) mcReader->Restart();

            // Declare TTreeReaderValues for data and MC for this layer
            TTreeReaderValue<double> cal_x(dataReader, x_branch.c_str());
            TTreeReaderValue<double> cal_y(dataReader, y_branch.c_str());
            TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

            TTreeReaderValue<double>* mc_cal_x = nullptr;
            TTreeReaderValue<double>* mc_cal_y = nullptr;
            TTreeReaderValue<int>* mc_particle_pid = nullptr;

            if (mcReader) {
                mc_cal_x = new TTreeReaderValue<double>(*mcReader, x_branch.c_str());
                mc_cal_y = new TTreeReaderValue<double>(*mcReader, y_branch.c_str());
                mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
            }

            // Create histograms for data and MC for each strictness level
            TH2D* h_data_0 = new TH2D("h_data_0", ("data " + layer_name + " hit position (" + particle_name + ", strictness_0)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
            TH2D* h_data_1 = new TH2D("h_data_1", ("data " + layer_name + " hit position (" + particle_name + ", strictness_1)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
            TH2D* h_data_2 = new TH2D("h_data_2", ("data " + layer_name + " hit position (" + particle_name + ", strictness_2)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
            TH2D* h_data_3 = new TH2D("h_data_3", ("data " + layer_name + " hit position (" + particle_name + ", strictness_3)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);

            h_data_0->GetXaxis()->SetTitle(("x_{" + layer_name + "}").c_str());
            h_data_0->GetYaxis()->SetTitle(("y_{" + layer_name + "}").c_str());
            h_data_1->GetXaxis()->SetTitle(("x_{" + layer_name + "}").c_str());
            h_data_1->GetYaxis()->SetTitle(("y_{" + layer_name + "}").c_str());
            h_data_2->GetXaxis()->SetTitle(("x_{" + layer_name + "}").c_str());
            h_data_2->GetYaxis()->SetTitle(("y_{" + layer_name + "}").c_str());
            h_data_3->GetXaxis()->SetTitle(("x_{" + layer_name + "}").c_str());
            h_data_3->GetYaxis()->SetTitle(("y_{" + layer_name + "}").c_str());

            TH2D* h_mc_0 = nullptr;
            TH2D* h_mc_1 = nullptr;
            TH2D* h_mc_2 = nullptr;
            TH2D* h_mc_3 = nullptr;

            if (mcReader) {
                h_mc_0 = new TH2D("h_mc_0", ("mc " + layer_name + " hit position (" + particle_name + ", strictness_0)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
                h_mc_1 = new TH2D("h_mc_1", ("mc " + layer_name + " hit position (" + particle_name + ", strictness_1)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
                h_mc_2 = new TH2D("h_mc_2", ("mc " + layer_name + " hit position (" + particle_name + ", strictness_2)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
                h_mc_3 = new TH2D("h_mc_3", ("mc " + layer_name + " hit position (" + particle_name + ", strictness_3)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);

                h_mc_0->GetXaxis()->SetTitle(("x_{" + layer_name + "}").c_str());
                h_mc_0->GetYaxis()->SetTitle(("y_{" + layer_name + "}").c_str());
                h_mc_1->GetXaxis()->SetTitle(("x_{" + layer_name + "}").c_str());
                h_mc_1->GetYaxis()->SetTitle(("y_{" + layer_name + "}").c_str());
                h_mc_2->GetXaxis()->SetTitle(("x_{" + layer_name + "}").c_str());
                h_mc_2->GetYaxis()->SetTitle(("y_{" + layer_name + "}").c_str());
                h_mc_3->GetXaxis()->SetTitle(("x_{" + layer_name + "}").c_str());
                h_mc_3->GetYaxis()->SetTitle(("y_{" + layer_name + "}").c_str());
            }

            // Fill the data histograms, applying the cuts
            while (dataReader.Next()) {
                if (*particle_pid == pid && *cal_x != -9999 && *cal_y != -9999) {
                    h_data_0->Fill(*cal_x, *cal_y); // No cuts
                    if (pcal_fiducial(*cal_lv_1, *cal_lw_1, *cal_lu_1, 
                    	*cal_lv_4, *cal_lw_4, *cal_lu_4, 
                    	*cal_lv_7, *cal_lw_7, *cal_lu_7, *cal_sector, 1)) {
                        h_data_1->Fill(*cal_x, *cal_y);
                    }
                    if (pcal_fiducial(*cal_lv_1, *cal_lw_1, *cal_lu_1, 
                    	*cal_lv_4, *cal_lw_4, *cal_lu_4, 
                    	*cal_lv_7, *cal_lw_7, *cal_lu_7, *cal_sector, 2)) {
                        h_data_2->Fill(*cal_x, *cal_y);
                    }
                    if (pcal_fiducial(*cal_lv_1, *cal_lw_1, *cal_lu_1, 
                    	*cal_lv_4, *cal_lw_4, *cal_lu_4, 
                    	*cal_lv_7, *cal_lw_7, *cal_lu_7, *cal_sector, 3)) {
                        h_data_3->Fill(*cal_x, *cal_y);
                    }
                }
            }

            // Fill the MC histograms if available, applying the cuts
            if (mcReader) {
                while (mcReader->Next()) {
                if (**mc_particle_pid == pid && **mc_cal_x != -9999 && **mc_cal_y != -9999) {
                    h_mc_0->Fill(**mc_cal_x, **mc_cal_y); // No cuts
                if (pcal_fiducial(**mc_cal_lv_1, **mc_cal_lw_1, **mc_cal_lu_1, 
                	**mc_cal_lv_4, **mc_cal_lw_4, **mc_cal_lu_4, 
                	**mc_cal_lv_7, **mc_cal_lw_7, **mc_cal_lu_7, 
                	**mc_cal_sector, 1)) {
                    h_mc_1->Fill(**mc_cal_x, **mc_cal_y);
                }
                if (pcal_fiducial(**mc_cal_lv_1, **mc_cal_lw_1, **mc_cal_lu_1, 
                	**mc_cal_lv_4, **mc_cal_lw_4, **mc_cal_lu_4, 
                	**mc_cal_lv_7, **mc_cal_lw_7, **mc_cal_lu_7, 
                	**mc_cal_sector, 2)) {
                    h_mc_2->Fill(**mc_cal_x, **mc_cal_y);
                }
                if (pcal_fiducial(**mc_cal_lv_1, **mc_cal_lw_1, **mc_cal_lu_1, 
                	**mc_cal_lv_4, **mc_cal_lw_4, **mc_cal_lu_4, 
                	**mc_cal_lv_7, **mc_cal_lw_7, **mc_cal_lu_7, 
                	**mc_cal_sector, 3)) {
                    h_mc_3->Fill(**mc_cal_x, **mc_cal_y);
                }
            }
        }
    }
    // Create a canvas to hold the 2x4 subplots
    TCanvas* c = new TCanvas(("c_" + particle_name + "_" + layer_name).c_str(), ("c_" + particle_name + "_" + layer_name).c_str(), 1600, 800);
    c->Divide(4, 2);

    // Draw the data plots on the top row
    for (int i = 1; i <= 4; ++i) {
        c->cd(i);
        gPad->SetLogz();             // Set log scale for the z-axis
        gPad->SetMargin(0.15, 0.15, 0.1, 0.1); // Increase padding
        if (i == 1) h_data_0->Draw("COLZ");
        if (i == 2) h_data_1->Draw("COLZ");
        if (i == 3) h_data_2->Draw("COLZ");
        if (i == 4) h_data_3->Draw("COLZ");
    }

    // Draw the MC plots on the bottom row, if available
    if (mcReader) {
        for (int i = 5; i <= 8; ++i) {
            c->cd(i);
            gPad->SetLogz();             // Set log scale for the z-axis
            gPad->SetMargin(0.15, 0.15, 0.1, 0.1); // Increase padding
            if (i == 5) h_mc_0->Draw("COLZ");
            if (i == 6) h_mc_1->Draw("COLZ");
            if (i == 7) h_mc_2->Draw("COLZ");
            if (i == 8) h_mc_3->Draw("COLZ");
        }
    }

    // Save the canvas
    c->SaveAs(("output/calibration/cal/" + particle_name + "_" + layer_name + "_cal_hit_position.png").c_str());

    // Clean up for this layer and particle type
    delete h_data_0;
    delete h_data_1;
    delete h_data_2;
    delete h_data_3;
    if (h_mc_0) delete h_mc_0;
    if (h_mc_1) delete h_mc_1;
            if (h_mc_2) delete h_mc_2;
            if (h_mc_3) delete h_mc_3;
            delete c;
            if (mc_cal_x) delete mc_cal_x;
            if (mc_cal_y) delete mc_cal_y;
            if (mc_particle_pid) delete mc_particle_pid;
        }
    }

    // Clean up the dynamically allocated memory for layer 1 variables
    if (mc_cal_lv_1) delete mc_cal_lv_1;
    if (mc_cal_lw_1) delete mc_cal_lw_1;
    if (mc_cal_lu_1) delete mc_cal_lu_1;
    if (mc_cal_sector) delete mc_cal_sector;
}

bool dc_fiducial(double edge_6, double edge_18, double edge_36, 
	int pid) {
    if (pid == 11 || pid == -11) {
        return edge_6 > 5 && edge_18 > 5 && edge_36 > 10;
    } else if (pid == 211 || pid == -211 || pid == 321 || pid == -321 || pid == 2212 || pid == -2212) {
        return edge_6 > 3 && edge_18 > 3 && edge_36 > 7;
    } 
    return false; // not a charged track? wrong pid?
}

void plot_dc_hit_position(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    int nBins = 100;

    std::vector<std::tuple<std::string, std::string, std::string, double, double>> regions = {
        {"traj_x_6", "traj_y_6", "region_1", -200, 200},
        {"traj_x_18", "traj_y_18", "region_2", -300, 300},
        {"traj_x_36", "traj_y_36", "region_3", -450, 450}
    };

    std::vector<std::tuple<int, std::string>> particle_types = {
        {11, "electron"},
        {-211, "pim"},
        {211, "pip"},
        {321, "kp"},
        {-321, "km"},
        {2212, "proton"}
    };

    // Declare TTreeReaderValues for the DC edge and track variables
    TTreeReaderValue<double> traj_edge_6(dataReader, "traj_edge_6");
    TTreeReaderValue<double> traj_edge_18(dataReader, "traj_edge_18");
    TTreeReaderValue<double> traj_edge_36(dataReader, "traj_edge_36");

    TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

    TTreeReaderValue<double>* mc_traj_edge_6 = nullptr;
    TTreeReaderValue<double>* mc_traj_edge_18 = nullptr;
    TTreeReaderValue<double>* mc_traj_edge_36 = nullptr;
    TTreeReaderValue<int>* mc_particle_pid = nullptr;

    if (mcReader) {
        mc_traj_edge_6 = new TTreeReaderValue<double>(*mcReader, "traj_edge_6");
        mc_traj_edge_18 = new TTreeReaderValue<double>(*mcReader, "traj_edge_18");
        mc_traj_edge_36 = new TTreeReaderValue<double>(*mcReader, "traj_edge_36");
        mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
    }

    // Declare TTreeReaderValues for trajectory x and y coordinates
    std::vector<TTreeReaderValue<double>> traj_x;
    std::vector<TTreeReaderValue<double>> traj_y;

    std::vector<TTreeReaderValue<double>*> mc_traj_x;
    std::vector<TTreeReaderValue<double>*> mc_traj_y;

    // Initialize TTreeReaderValues for each region
    for (const auto& region : regions) {
        traj_x.emplace_back(dataReader, std::get<0>(region).c_str());
        traj_y.emplace_back(dataReader, std::get<1>(region).c_str());

        if (mcReader) {
            mc_traj_x.push_back(new TTreeReaderValue<double>(*mcReader, std::get<0>(region).c_str()));
            mc_traj_y.push_back(new TTreeReaderValue<double>(*mcReader, std::get<1>(region).c_str()));
        }
    }

    for (const auto& particle_type : particle_types) {
        int pid = std::get<0>(particle_type);
        std::string particle_name = std::get<1>(particle_type);

        // Create a canvas for data
        TCanvas* c_data = new TCanvas(("c_data_" + particle_name).c_str(), ("Data DC Hit Position (" + particle_name + ")").c_str(), 1800, 1200);
        c_data->Divide(3, 2);

        TCanvas* c_mc = nullptr;
        if (mcReader) {
            c_mc = new TCanvas(("c_mc_" + particle_name).c_str(), ("MC DC Hit Position (" + particle_name + ")").c_str(), 1800, 1200);
            c_mc->Divide(3, 2);
        }

        // Create histograms for data and MC
        std::vector<TH2D*> h_data_before(3), h_data_after(3);
        std::vector<TH2D*> h_mc_before(3), h_mc_after(3);

        for (int region_idx = 0; region_idx < 3; ++region_idx) {
            std::string region_name = std::get<2>(regions[region_idx]);
            double xMin = std::get<3>(regions[region_idx]);
            double xMax = std::get<4>(regions[region_idx]);
            double yMin = xMin;
            double yMax = xMax;

            h_data_before[region_idx] = new TH2D(("h_data_before_" + region_name).c_str(), ("Data " + region_name + " Before Cuts (" + particle_name + ")").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
            h_data_after[region_idx] = new TH2D(("h_data_after_" + region_name).c_str(), ("Data " + region_name + " After Cuts (" + particle_name + ")").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);

            if (mcReader) {
                h_mc_before[region_idx] = new TH2D(("h_mc_before_" + region_name).c_str(), ("MC " + region_name + " Before Cuts (" + particle_name + ")").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
                h_mc_after[region_idx] = new TH2D(("h_mc_after_" + region_name).c_str(), ("MC " + region_name + " After Cuts (" + particle_name + ")").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
            }
        }

        // Fill the data histograms
        dataReader.Restart();
        while (dataReader.Next()) {
            if (*particle_pid == pid) {
                for (int region_idx = 0; region_idx < 3; ++region_idx) {
                    double traj_x_value = *traj_x[region_idx];
                    double traj_y_value = *traj_y[region_idx];

                    if (traj_x_value != -9999 && traj_y_value != -9999) {
                        h_data_before[region_idx]->Fill(traj_x_value, traj_y_value);
                        if (dc_fiducial(*traj_edge_6, *traj_edge_18, *traj_edge_36, pid)) {
                            h_data_after[region_idx]->Fill(traj_x_value, traj_y_value);
                        }
                    }
                }
            }
        }

        // Fill the MC histograms if available
        if (mcReader) {
            mcReader->Restart();
            while (mcReader->Next()) {
                if (**mc_particle_pid == pid) {
                    for (int region_idx = 0; region_idx < 3; ++region_idx) {
                        double mc_traj_x_value = **mc_traj_x[region_idx];
                        double mc_traj_y_value = **mc_traj_y[region_idx];

                        if (mc_traj_x_value != -9999 && mc_traj_y_value != -9999) {
                            h_mc_before[region_idx]->Fill(mc_traj_x_value, mc_traj_y_value);
                            if (dc_fiducial(**mc_traj_edge_6, **mc_traj_edge_18, **mc_traj_edge_36, pid)) {
                                h_mc_after[region_idx]->Fill(mc_traj_x_value, mc_traj_y_value);
                            }
                        }
                    }
                }
            }
        }

        // Find the maximum value across all histograms for consistent scaling
        double max_value_data = 0, max_value_mc = 0;
        for (int region_idx = 0; region_idx < 3; ++region_idx) {
            max_value_data = std::max(max_value_data, h_data_before[region_idx]->GetMaximum());
            max_value_data = std::max(max_value_data, h_data_after[region_idx]->GetMaximum());
            if (mcReader) {
                max_value_mc = std::max(max_value_mc, h_mc_before[region_idx]->GetMaximum());
                max_value_mc = std::max(max_value_mc, h_mc_after[region_idx]->GetMaximum());
            }
        }

        // Set the maximum for each histogram to ensure consistent scaling
        for (int region_idx = 0; region_idx < 3; ++region_idx) {
            h_data_before[region_idx]->SetMaximum(max_value_data * 1.1);
            h_data_after[region_idx]->SetMaximum(max_value_data * 1.1);
            if (mcReader) {
                h_mc_before[region_idx]->SetMaximum(max_value_mc * 1.1);
                h_mc_after[region_idx]->SetMaximum(max_value_mc * 1.1);
            }
        }

        // Draw and save the data canvas
        for (int region_idx = 0; region_idx < 3; ++region_idx) {
            c_data->cd(region_idx + 1);
            gPad->SetLogz();
            gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
            h_data_before[region_idx]->Draw("COLZ");
            c_data->cd(region_idx + 4);
            gPad->SetLogz();
            gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
            h_data_after[region_idx]->Draw("COLZ");
        }
        c_data->SaveAs(("output/calibration/dc/positions/data_" + particle_name + "_dc_hit_position.png").c_str());

        // Draw and save the MC canvas if available
        if (mcReader) {
            for (int region_idx = 0; region_idx < 3; ++region_idx) {
                c_mc->cd(region_idx + 1);
                gPad->SetLogz();
                gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
                h_mc_before[region_idx]->Draw("COLZ");
                c_mc->cd(region_idx + 4);
                gPad->SetLogz();
                gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
                h_mc_after[region_idx]->Draw("COLZ");
            }
            c_mc->SaveAs(("output/calibration/dc/positions/mc_" + particle_name + "_dc_hit_position.png").c_str());
        }

        // Cleanup
        for (int region_idx = 0; region_idx < 3; ++region_idx) {
            delete h_data_before[region_idx];
            delete h_data_after[region_idx];
            if (mcReader) {
                delete h_mc_before[region_idx];
                delete h_mc_after[region_idx];
            }
        }
        delete c_data;
        if (mcReader) delete c_mc;
    }

    // Clean up the dynamically allocated memory for edge variables
    if (mc_traj_edge_6) delete mc_traj_edge_6;
    if (mc_traj_edge_18) delete mc_traj_edge_18;
    if (mc_traj_edge_36) delete mc_traj_edge_36;

    // Clean up mc_traj_x and mc_traj_y vectors
    for (auto& ptr : mc_traj_x) delete ptr;
    for (auto& ptr : mc_traj_y) delete ptr;
}

void normalize_histogram(TH2D* h_sum, TH2D* h_count) {
    for (int i = 1; i <= h_sum->GetNbinsX(); ++i) {
        for (int j = 1; j <= h_sum->GetNbinsY(); ++j) {
            double count = h_count->GetBinContent(i, j);
            if (count > 0) {
                h_sum->SetBinContent(i, j, h_sum->GetBinContent(i, j) / count);
            } else {
                h_sum->SetBinContent(i, j, 0);
            }
        }
    }
}

std::vector<TH2D*> create_histograms_for_sector(const std::string& region_name, const std::string& particle_name, int nBins, double xMin, double xMax, double yMin, double yMax, bool isMC) {
    std::vector<TH2D*> histograms(6);

    for (int sector = 0; sector < 6; ++sector) {
        histograms[sector] = new TH2D(
            ((isMC ? "h_mc_sum_" : "h_data_sum_") + region_name + "_sector" + std::to_string(sector + 1)).c_str(),
            ((isMC ? "mc " : "data ") + region_name + " sector " + std::to_string(sector + 1) + " #chi^{2}/ndf (" + particle_name + ")").c_str(),
            nBins, xMin, xMax, nBins, yMin, yMax
        );
        histograms[sector]->SetDirectory(0);  // Detach histogram from current ROOT directory
    }

    return histograms;
}

void draw_and_save_sector_histograms(TCanvas* canvas, std::vector<TH2D*>& histograms, const std::string& output_file) {
    double max_value = 0;

    // Find the maximum z-axis value across all histograms
    for (const auto& hist : histograms) {
        double hist_max = hist->GetMaximum();
        if (hist_max > max_value) {
            max_value = hist_max;
        }
    }

    // Set the maximum for each histogram and draw
    for (int sector = 0; sector < 6; ++sector) {
        canvas->cd(sector + 1);  // Select the pad corresponding to the sector
        gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
        gPad->SetLogz();  // Set log scale for the z-axis
        histograms[sector]->SetStats(false);  // Disable stat box
        histograms[sector]->SetMaximum(max_value);  // Set the same max value for z-axis
        histograms[sector]->Draw("COLZ");
    }
    canvas->SaveAs(output_file.c_str());
}

std::pair<double, double> rotate_coordinates(double x, double y, int sector) {
    double angle = -60.0 * (sector - 1);  // Angle to rotate counterclockwise for sectors 2-6
    double radians = angle * TMath::Pi() / 180.0;

    double x_rot = x * cos(radians) - y * sin(radians);
    double y_rot = x * sin(radians) + y * cos(radians);

    return std::make_pair(x_rot, y_rot);
}

void dc_fiducial_determination(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    int nBins = 100;
    std::vector<std::tuple<std::string, std::string, std::string, double, double, double, double, std::string, double>> regions = {
        {"traj_x_6", "traj_y_6", "region_1", 15, 160, -80, 80, "traj_edge_6", 30},
        {"traj_x_18", "traj_y_18", "region_2", 30, 240, -125, 125, "traj_edge_18", 30},
        {"traj_x_36", "traj_y_36", "region_3", 30, 400, -200, 200, "traj_edge_36", 30}
    };

    // Array of particle types (photons and electrons) and their corresponding PIDs
    std::vector<std::tuple<int, std::string>> particle_types = {
        // {11, "electron"},
        // {-211, "pim"},
        // {211, "pip"},
        // {321, "kp"},
        // {-321, "km"},
        {2212, "proton"}
    };

    for (const auto& particle_type : particle_types) {
        int pid = std::get<0>(particle_type);
        std::string particle_name = std::get<1>(particle_type);

        for (const auto& region : regions) {
            std::string x_branch = std::get<0>(region);
            std::string y_branch = std::get<1>(region);
            std::string region_name = std::get<2>(region);
            double xMin = std::get<3>(region);
            double xMax = std::get<4>(region);
            double yMin = std::get<5>(region);
            double yMax = std::get<6>(region);
            std::string edge_branch = std::get<7>(region);
            double edge_max = std::get<8>(region);

            // Create histograms for data and MC
            auto h_data_sum_sector = create_histograms_for_sector(region_name, particle_name, nBins, xMin, xMax, yMin, yMax, false);
            auto h_data_count_sector = create_histograms_for_sector(region_name, particle_name, nBins, xMin, xMax, yMin, yMax, false);

            auto h_mc_sum_sector = create_histograms_for_sector(region_name, particle_name, nBins, xMin, xMax, yMin, yMax, true);
            auto h_mc_count_sector = create_histograms_for_sector(region_name, particle_name, nBins, xMin, xMax, yMin, yMax, true);

            TCanvas* c_region = new TCanvas(("c_" + particle_name + "_" + region_name + "_chi2_ndf").c_str(), ("c_" + particle_name + " #chi^{2}/ndf").c_str(), 1800, 1200);
            c_region->Divide(3, 2);
            TCanvas* c_mc_region = nullptr;
            if (mcReader) {
                c_mc_region = new TCanvas(("c_mc_" + particle_name + "_" + region_name + "_chi2_ndf").c_str(), ("MC " + particle_name + " #chi^{2}/ndf").c_str(), 1800, 1200);
                c_mc_region->Divide(3, 2);
            }

            // Restart readers
            dataReader.Restart();
            if (mcReader) mcReader->Restart();

            TTreeReaderValue<double> traj_x(dataReader, x_branch.c_str());
            TTreeReaderValue<double> traj_y(dataReader, y_branch.c_str());
            TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");
            TTreeReaderValue<int> track_sector_6(dataReader, "track_sector_6");
            TTreeReaderValue<double> track_chi2_6(dataReader, "track_chi2_6");
            TTreeReaderValue<int> track_ndf_6(dataReader, "track_ndf_6");
            TTreeReaderValue<double> traj_edge(dataReader, edge_branch.c_str());
            TTreeReaderValue<double> traj_edge_6(dataReader, "traj_edge_6");
            TTreeReaderValue<double> traj_edge_18(dataReader, "traj_edge_18");
            TTreeReaderValue<double> traj_edge_36(dataReader, "traj_edge_36");
            TTreeReaderValue<double> track_theta(dataReader, "theta");

            TTreeReaderValue<double>* mc_traj_x = nullptr;
            TTreeReaderValue<double>* mc_traj_y = nullptr;
            TTreeReaderValue<int>* mc_particle_pid = nullptr;
            TTreeReaderValue<int>* mc_track_sector_6 = nullptr;
            TTreeReaderValue<double>* mc_track_chi2_6 = nullptr;
            TTreeReaderValue<int>* mc_track_ndf_6 = nullptr;
            TTreeReaderValue<double>* mc_traj_edge = nullptr;
            TTreeReaderValue<double>* mc_traj_edge_6 = nullptr;
            TTreeReaderValue<double>* mc_traj_edge_18 = nullptr;
            TTreeReaderValue<double>* mc_traj_edge_36 = nullptr;
            TTreeReaderValue<double>* mc_track_theta = nullptr;

            if (mcReader) {
                mc_traj_x = new TTreeReaderValue<double>(*mcReader, x_branch.c_str());
                mc_traj_y = new TTreeReaderValue<double>(*mcReader, y_branch.c_str());
                mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
                mc_track_sector_6 = new TTreeReaderValue<int>(*mcReader, "track_sector_6");
                mc_track_chi2_6 = new TTreeReaderValue<double>(*mcReader, "track_chi2_6");
                mc_track_ndf_6 = new TTreeReaderValue<int>(*mcReader, "track_ndf_6");
                mc_traj_edge = new TTreeReaderValue<double>(*mcReader, edge_branch.c_str());
                mc_traj_edge_6 = new TTreeReaderValue<double>(*mcReader, "traj_edge_6");
                mc_traj_edge_18 = new TTreeReaderValue<double>(*mcReader, "traj_edge_18");
                mc_traj_edge_36 = new TTreeReaderValue<double>(*mcReader, "traj_edge_36");
                mc_track_theta = new TTreeReaderValue<double>(*mcReader, "theta");
            }

            // Fill data histograms for the rotated hit positions
            while (dataReader.Next()) {
                if (*particle_pid == pid && *traj_x != -9999 && *traj_y != -9999 && *track_ndf_6 > 0) {
                    double chi2_ndf = *track_chi2_6 / *track_ndf_6;

                    auto rotated_coords = rotate_coordinates(*traj_x, *traj_y, *track_sector_6);
                    double traj_x_rot = rotated_coords.first;
                    double traj_y_rot = rotated_coords.second;

                    if (traj_x_rot >= 0) {
                        h_data_sum_sector[*track_sector_6 - 1]->Fill(traj_x_rot, traj_y_rot, chi2_ndf);
                        h_data_count_sector[*track_sector_6 - 1]->Fill(traj_x_rot, traj_y_rot);
                    }
                }
            }

            if (mcReader) {
                while (mcReader->Next()) {
                    if (**mc_particle_pid == pid && **mc_traj_x != -9999 && **mc_traj_y != -9999 && **mc_track_ndf_6 > 0) {
                        double mc_chi2_ndf = **mc_track_chi2_6 / **mc_track_ndf_6;

                        auto rotated_coords = rotate_coordinates(**mc_traj_x, **mc_traj_y, **mc_track_sector_6);
                        double mc_traj_x_rot = rotated_coords.first;
                        double mc_traj_y_rot = rotated_coords.second;

                        if (mc_traj_x_rot >= 0) {
                            h_mc_sum_sector[**mc_track_sector_6 - 1]->Fill(mc_traj_x_rot, mc_traj_y_rot, mc_chi2_ndf);
                            h_mc_count_sector[**mc_track_sector_6 - 1]->Fill(mc_traj_x_rot, mc_traj_y_rot);
                        }
                    }
                }
            }

            // Normalize and save rotated hit position histograms
            for (int sector = 0; sector < 6; ++sector) {
                normalize_histogram(h_data_sum_sector[sector], h_data_count_sector[sector]);

                if (mcReader) {
                    normalize_histogram(h_mc_sum_sector[sector], h_mc_count_sector[sector]);
                }
            }

            draw_and_save_sector_histograms(c_region, h_data_sum_sector, "output/calibration/dc/determination/chi2_per_ndf_" + particle_name + "_" + region_name + ".png");

            if (mcReader) {
                draw_and_save_sector_histograms(c_mc_region, h_mc_sum_sector, "output/calibration/dc/determination/mc_chi2_per_ndf_" + particle_name + "_" + region_name + ".png");
            }

            // Cleanup for the rotated hit position histograms
            delete c_region;
            if (c_mc_region) delete c_mc_region;

            for (auto& hist : h_data_sum_sector) delete hist;
            for (auto& hist : h_data_count_sector) delete hist;

            if (mcReader) {
                for (auto& hist : h_mc_sum_sector) delete hist;
                for (auto& hist : h_mc_count_sector) delete hist;
            }

            // Now, calculate and plot mean chi2/ndf as a function of traj_edge for each sector
            std::vector<TH1D*> h_sum_chi2_ndf_sector(6);
            std::vector<TH1D*> h_count_chi2_ndf_sector(6);
            std::vector<TH1D*> h_sum_chi2_ndf_mc_sector(6);
            std::vector<TH1D*> h_count_chi2_ndf_mc_sector(6);
            for (int sector = 0; sector < 6; ++sector) {
                h_sum_chi2_ndf_sector[sector] = new TH1D(("h_sum_chi2_ndf_sector_" + region_name + "_sector" + std::to_string(sector + 1)).c_str(),
                                                         (particle_name + " in " + region_name + " - Sector " + std::to_string(sector + 1)).c_str(),
                                                         nBins, 0, edge_max);
                h_sum_chi2_ndf_sector[sector]->GetXaxis()->SetTitle("edge");
                h_sum_chi2_ndf_sector[sector]->GetYaxis()->SetTitle("<chi2/ndf>");
                h_count_chi2_ndf_sector[sector] = new TH1D(("h_count_chi2_ndf_sector_" + region_name + "_sector" + std::to_string(sector + 1)).c_str(),
                                                       "", nBins, 0, edge_max);

                if (mcReader) {
                    h_sum_chi2_ndf_mc_sector[sector] = new TH1D(("h_sum_chi2_ndf_mc_sector_" + region_name + "_sector" + std::to_string(sector + 1)).c_str(),
                                                                (particle_name + " in " + region_name + " - MC - Sector " + std::to_string(sector + 1)).c_str(),
                                                                nBins, 0, edge_max);
                    h_sum_chi2_ndf_mc_sector[sector]->GetXaxis()->SetTitle("edge");
                    h_sum_chi2_ndf_mc_sector[sector]->GetYaxis()->SetTitle("<chi2/ndf>");

                    h_count_chi2_ndf_mc_sector[sector] = new TH1D(("h_count_chi2_ndf_mc_sector_" + region_name + "_sector" + std::to_string(sector + 1)).c_str(),
                                                                  "", nBins, 0, edge_max);
                }
            }

            // Fill histograms for data
            dataReader.Restart();
            while (dataReader.Next()) {
                if (*particle_pid == pid && *traj_edge != -9999 && *track_ndf_6 > 0) {
                    double chi2_ndf = *track_chi2_6 / *track_ndf_6;
                    h_sum_chi2_ndf_sector[*track_sector_6 - 1]->Fill(*traj_edge, chi2_ndf);
                    h_count_chi2_ndf_sector[*track_sector_6 - 1]->Fill(*traj_edge);
                }
            }

            // Fill histograms for MC
            if (mcReader) {
                mcReader->Restart();
                while (mcReader->Next()) {
                    if (**mc_particle_pid == pid && **mc_traj_edge != -9999 && **mc_track_ndf_6 > 0) {
                        double mc_chi2_ndf = **mc_track_chi2_6 / **mc_track_ndf_6;
                        h_sum_chi2_ndf_mc_sector[**mc_track_sector_6 - 1]->Fill(**mc_traj_edge, mc_chi2_ndf);
                        h_count_chi2_ndf_mc_sector[**mc_track_sector_6 - 1]->Fill(**mc_traj_edge);
                    }
                }
            }

            // Normalize the histograms to get the mean chi2/ndf per bin of traj_edge
            for (int sector = 0; sector < 6; ++sector) {
                h_sum_chi2_ndf_sector[sector]->Divide(h_count_chi2_ndf_sector[sector]);
                if (mcReader) {
                    h_sum_chi2_ndf_mc_sector[sector]->Divide(h_count_chi2_ndf_mc_sector[sector]);
                }
            }

            // Find the max y-value for setting axis limits
            double max_value = 0;
            for (int sector = 0; sector < 6; ++sector) {
                if (h_sum_chi2_ndf_sector[sector]->GetMaximum() > max_value) {
                    max_value = h_sum_chi2_ndf_sector[sector]->GetMaximum();
                }
                if (mcReader && h_sum_chi2_ndf_mc_sector[sector]->GetMaximum() > max_value) {
                    max_value = h_sum_chi2_ndf_mc_sector[sector]->GetMaximum();
                }
            }
            max_value *= 1.1; // Add 10% margin to max value

            // Draw and save the mean chi2/ndf vs traj_edge plots
            TCanvas* c_edge = new TCanvas(("c_edge_" + particle_name + "_" + region_name).c_str(), ("Mean chi2/ndf vs traj_edge for " + particle_name + " in " + region_name).c_str(), 1800, 1200);
            c_edge->Divide(3, 2);

            for (int sector = 0; sector < 6; ++sector) {
                c_edge->cd(sector + 1);
                gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
                gPad->SetLogy();  
                h_sum_chi2_ndf_sector[sector]->SetStats(false);
                h_sum_chi2_ndf_sector[sector]->SetAxisRange(0, 30, "X"); 
                // h_sum_chi2_ndf_sector[sector]->SetMaximum(max_value);
                h_sum_chi2_ndf_sector[sector]->SetMaximum(1e3);
                h_sum_chi2_ndf_sector[sector]->SetMinimum(10);
                h_sum_chi2_ndf_sector[sector]->SetLineColor(kBlack);
                h_sum_chi2_ndf_sector[sector]->Draw("E");

                if (mcReader) {
                    h_sum_chi2_ndf_mc_sector[sector]->SetStats(false);
                    h_sum_chi2_ndf_mc_sector[sector]->SetAxisRange(0, 30, "X"); 
                    // h_sum_chi2_ndf_mc_sector[sector]->SetMaximum(max_value);
                    h_sum_chi2_ndf_mc_sector[sector]->SetMaximum(1e3);
                    h_sum_chi2_ndf_mc_sector[sector]->SetMinimum(10);
                    h_sum_chi2_ndf_mc_sector[sector]->SetLineColor(kRed);
                    h_sum_chi2_ndf_mc_sector[sector]->Draw("E SAME");
                }

                TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
                legend->AddEntry(h_sum_chi2_ndf_sector[sector], "Data", "l");
                if (mcReader) legend->AddEntry(h_sum_chi2_ndf_mc_sector[sector], "MC", "l");
                legend->Draw();
            }

            c_edge->SaveAs(("output/calibration/dc/determination/mean_chi2_per_ndf_vs_traj_edge_" + particle_name + "_" + region_name + ".png").c_str());

            // Cleanup for mean chi2/ndf vs traj_edge histograms
            delete c_edge;
            for (auto& hist : h_sum_chi2_ndf_sector) delete hist;
            for (auto& hist : h_count_chi2_ndf_sector) delete hist;

            if (mcReader) {
                for (auto& hist : h_sum_chi2_ndf_mc_sector) delete hist;
                for (auto& hist : h_count_chi2_ndf_mc_sector) delete hist;
            }

            if (mc_traj_edge) delete mc_traj_edge;

            // Count the number of non-kaon particles
            int num_valid_particles = 0;
            int num_particles = particle_types.size();
            for (int particle_idx = 0; particle_idx < num_particles; ++particle_idx) {
                int pid = std::get<0>(particle_types[particle_idx]);
                if (pid != 321 && pid != -321) {
                    ++num_valid_particles;
                }
            }

            // Create canvas with correct dimensions based on non-kaon particles
            TCanvas* c_theta = new TCanvas("c_theta", "Mean chi2/ndf vs Theta", 1800 * num_valid_particles, 1200);
            c_theta->Divide(num_valid_particles, 2);

            std::vector<TH2D*> h2_chi2_vs_theta_data(num_valid_particles);
            std::vector<TH2D*> h2_chi2_vs_theta_mc(num_valid_particles);

            int valid_particle_idx = 0;
            for (int particle_idx = 0; particle_idx < num_particles; ++particle_idx) {
                int pid = std::get<0>(particle_types[particle_idx]);

                // Skip kaons
                if (pid == 321 || pid == -321) {
                    continue;
                }

                std::string particle_name = std::get<1>(particle_types[particle_idx]);

                h2_chi2_vs_theta_data[valid_particle_idx] = new TH2D(
                    ("h2_chi2_vs_theta_data_" + particle_name).c_str(),
                    (particle_name + " (Data)").c_str(),
                    nBins, 0, 50, nBins, 0, 100
                );

                if (mcReader) {
                    h2_chi2_vs_theta_mc[valid_particle_idx] = new TH2D(
                        ("h2_chi2_vs_theta_mc_" + particle_name).c_str(),
                        (particle_name + " (MC)").c_str(),
                        nBins, 0, 50, nBins, 0, 100
                    );
                }

                // Fill the histograms with fiducial cuts applied for data
                dataReader.Restart();
                while (dataReader.Next()) {
                    if (*particle_pid == pid && *track_ndf_6 > 0 && 
                        dc_fiducial(*traj_edge_6, *traj_edge_18, *traj_edge_36, pid)) {
                        double chi2_ndf = *track_chi2_6 / *track_ndf_6;
                        h2_chi2_vs_theta_data[valid_particle_idx]->Fill(*track_theta, chi2_ndf);
                    }
                }

                // Fill the histograms with fiducial cuts applied for MC
                if (mcReader) {
                    mcReader->Restart();
                    while (mcReader->Next()) {
                        if (**mc_particle_pid == pid && **mc_track_ndf_6 > 0 && 
                            dc_fiducial(**mc_traj_edge_6, **mc_traj_edge_18, **mc_traj_edge_36, pid)) {
                            double mc_chi2_ndf = **mc_track_chi2_6 / **mc_track_ndf_6;
                            h2_chi2_vs_theta_mc[valid_particle_idx]->Fill(**mc_track_theta, mc_chi2_ndf);
                        }
                    }
                }

                // Draw and save the 2D histograms of chi2/ndf vs theta
                c_theta->cd(valid_particle_idx + 1);
                gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
                gPad->SetRightMargin(0.2);  // Adjust the right margin
                gPad->SetLogz();  // Set log scale for the z-axis
                h2_chi2_vs_theta_data[valid_particle_idx]->SetStats(false);
                h2_chi2_vs_theta_data[valid_particle_idx]->GetXaxis()->SetTitle("#theta (degrees)");
                h2_chi2_vs_theta_data[valid_particle_idx]->GetYaxis()->SetTitle("<chi2/ndf>");
                h2_chi2_vs_theta_data[valid_particle_idx]->Draw("COLZ");

                if (mcReader) {
                    c_theta->cd(valid_particle_idx + 1 + num_valid_particles); // Second row for MC
                    gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
                    gPad->SetRightMargin(0.2);  // Adjust the right margin
                    gPad->SetLogz();  // Set log scale for the z-axis
                    h2_chi2_vs_theta_mc[valid_particle_idx]->SetStats(false);
                    h2_chi2_vs_theta_mc[valid_particle_idx]->GetXaxis()->SetTitle("#theta (degrees)");
                    h2_chi2_vs_theta_mc[valid_particle_idx]->GetYaxis()->SetTitle("<chi2/ndf>");
                    h2_chi2_vs_theta_mc[valid_particle_idx]->Draw("COLZ");
                }

                ++valid_particle_idx;
            }

            // Save as PDF instead of PNG
            c_theta->SaveAs("output/calibration/dc/determination/mean_chi2_per_ndf_vs_theta_all_particles.png");

            // Cleanup
            for (int particle_idx = 0; particle_idx < num_valid_particles; ++particle_idx) {
                delete h2_chi2_vs_theta_data[particle_idx];
                if (mcReader) delete h2_chi2_vs_theta_mc[particle_idx];
            }
            delete c_theta;
        }

        dataReader.Restart();
        if (mcReader) mcReader->Restart();
    }
}

void plot_chi2pid_cd(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    int nBins = 100;
    double chi2pidMin = -10;
    double chi2pidMax = 10;
    double pMin = 0;
    double pMax = 8;  // Updated maximum momentum value
    double betaMax = 1.4;  // Updated maximum beta value
    std::vector<double> pBins = {0, 0.33, 0.67, 1.00, 1.33, 1.67, 2.00, 2.33, 2.67, 3.00, 3.5, 4, 4.5, 5, 6, 7, 8};

    // Particle types to analyze
    std::vector<std::tuple<int, std::string>> particle_types = {
        {211, "#pi^{+}"},
        {321, "k^{+}"},
        {2212, "p"},
        {-211, "#pi^{-}"},
        {-321, "k^{-}"},
        {-2212, "pbar"}
    };

    // Initialize readers for data and MC
    TTreeReaderValue<double> particle_chi2pid(dataReader, "particle_chi2pid");
    TTreeReaderValue<double> particle_p(dataReader, "p");
    TTreeReaderValue<double> particle_beta(dataReader, "particle_beta");  // Beta variable
    TTreeReaderValue<int> track_sector_6(dataReader, "track_sector_6");
    TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

    TTreeReaderValue<double>* mc_particle_chi2pid = nullptr;
    TTreeReaderValue<double>* mc_particle_p = nullptr;
    TTreeReaderValue<double>* mc_particle_beta = nullptr;  // MC Beta variable
    TTreeReaderValue<int>* mc_track_sector_6 = nullptr;
    TTreeReaderValue<int>* mc_particle_pid = nullptr;

    if (mcReader) {
        mc_particle_chi2pid = new TTreeReaderValue<double>(*mcReader, "particle_chi2pid");
        mc_particle_p = new TTreeReaderValue<double>(*mcReader, "p");
        mc_particle_beta = new TTreeReaderValue<double>(*mcReader, "particle_beta");
        mc_track_sector_6 = new TTreeReaderValue<int>(*mcReader, "track_sector_6");
        mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
    }

    // 1D Histograms canvas
    TCanvas* c = new TCanvas("c_chi2pid_cd", "chi2pid in CD", 1800, 1200);
    c->Divide(3, 2);
    gPad->SetLeftMargin(0.15);  // Add padding to the left

    // 2D Histograms canvas for data
    TCanvas* c_data_2D = new TCanvas("c_data_2D_chi2pid_vs_p_cd", "chi2pid vs p in CD", 1800, 1200);
    c_data_2D->Divide(3, 2);
    gPad->SetLeftMargin(0.15);  // Add padding to the left

    // 2D Histograms canvas for MC (if applicable)
    TCanvas* c_mc_2D = nullptr;
    if (mcReader) {
        c_mc_2D = new TCanvas("c_mc_2D_chi2pid_vs_p_cd", "chi2pid vs p in CD", 1800, 1200);
        c_mc_2D->Divide(3, 2);
        gPad->SetLeftMargin(0.15);  // Add padding to the left
    }

    // Additional 1x2 canvases for positive and negative tracks for particle_beta vs. p
    TCanvas* c_data_pos_neg_beta = new TCanvas("c_data_pos_neg_beta_vs_p_cd", "#beta vs p (Data)", 1800, 600);
    c_data_pos_neg_beta->Divide(2, 1);
    TCanvas* c_mc_pos_neg_beta = nullptr;
    if (mcReader) {
        c_mc_pos_neg_beta = new TCanvas("c_mc_pos_neg_beta_vs_p_cd", "#beta vs p (MC)", 1800, 600);
        c_mc_pos_neg_beta->Divide(2, 1);
    }

    // Initialize histograms for combined positive and negative tracks
    TH2D* h_data_beta_vs_p_pos = new TH2D("h_data_beta_vs_p_pos", "#beta vs p (Positive Tracks)", nBins, pMin, pMax, nBins, 0, betaMax);
    h_data_beta_vs_p_pos->GetXaxis()->SetTitle("p (GeV)");
    h_data_beta_vs_p_pos->GetYaxis()->SetTitle("#beta");
    h_data_beta_vs_p_pos->SetStats(false);

    TH2D* h_data_beta_vs_p_neg = new TH2D("h_data_beta_vs_p_neg", "#beta vs p (Negative Tracks)", nBins, pMin, pMax, nBins, 0, betaMax);
    h_data_beta_vs_p_neg->GetXaxis()->SetTitle("p (GeV)");
    h_data_beta_vs_p_neg->GetYaxis()->SetTitle("#beta");
    h_data_beta_vs_p_neg->SetStats(false);

    TH2D* h_mc_beta_vs_p_pos = nullptr;
    TH2D* h_mc_beta_vs_p_neg = nullptr;

    if (mcReader) {
        h_mc_beta_vs_p_pos = new TH2D("h_mc_beta_vs_p_pos", "#beta vs p (Positive Tracks)", nBins, pMin, pMax, nBins, 0, betaMax);
        h_mc_beta_vs_p_pos->GetXaxis()->SetTitle("p (GeV)");
        h_mc_beta_vs_p_pos->GetYaxis()->SetTitle("#beta");
        h_mc_beta_vs_p_pos->SetStats(false);

        h_mc_beta_vs_p_neg = new TH2D("h_mc_beta_vs_p_neg", "#beta vs p (Negative Tracks)", nBins, pMin, pMax, nBins, 0, betaMax);
        h_mc_beta_vs_p_neg->GetXaxis()->SetTitle("p (GeV)");
        h_mc_beta_vs_p_neg->GetYaxis()->SetTitle("#beta");
        h_mc_beta_vs_p_neg->SetStats(false);
    }

    // Initialize histograms for each particle type
    std::vector<TH1D*> h_data(6);
    std::vector<TH1D*> h_mc(6);
    std::vector<TH2D*> h_data_chi2pid_vs_p(6);
    std::vector<TH2D*> h_mc_chi2pid_vs_p(6);

    // Initialize histograms for beta distributions in momentum bins (4x4 canvases)
    std::vector<TH1D*> h_data_beta_bins_pos(pBins.size() - 1);
    std::vector<TH1D*> h_data_beta_bins_neg(pBins.size() - 1);
    std::vector<TH1D*> h_mc_beta_bins_pos(pBins.size() - 1);
    std::vector<TH1D*> h_mc_beta_bins_neg(pBins.size() - 1);

    for (size_t i = 0; i < particle_types.size(); ++i) {
        std::string hname = "h_data_" + std::get<1>(particle_types[i]);
        h_data[i] = new TH1D(hname.c_str(), (std::get<1>(particle_types[i])).c_str(), nBins, chi2pidMin, chi2pidMax);
        h_data[i]->GetXaxis()->SetTitle("chi2pid");
        h_data[i]->GetYaxis()->SetTitle("Normalized Counts");
        h_data[i]->SetStats(false);  // Hide the stat box

        if (mcReader) {
            hname = "h_mc_" + std::get<1>(particle_types[i]);
            h_mc[i] = new TH1D(hname.c_str(), (std::get<1>(particle_types[i])).c_str(), nBins, chi2pidMin, chi2pidMax);
            h_mc[i]->GetXaxis()->SetTitle("chi2pid");
            h_mc[i]->GetYaxis()->SetTitle("Normalized Counts");
            h_mc[i]->SetLineColor(kRed);
            h_mc[i]->SetStats(false);  // Hide the stat box
        }

        hname = "h_data_chi2pid_vs_p_" + std::get<1>(particle_types[i]);
        h_data_chi2pid_vs_p[i] = new TH2D(hname.c_str(), ("chi2pid vs p: " + std::get<1>(particle_types[i])).c_str(),
                                          nBins, pMin, pMax, nBins, chi2pidMin, chi2pidMax);
        h_data_chi2pid_vs_p[i]->GetXaxis()->SetTitle("p (GeV)");
        h_data_chi2pid_vs_p[i]->GetYaxis()->SetTitle("chi2pid");
        h_data_chi2pid_vs_p[i]->SetStats(false);  // Hide the stat box
        if (mcReader) {
            hname = "h_mc_chi2pid_vs_p_" + std::get<1>(particle_types[i]);
            h_mc_chi2pid_vs_p[i] = new TH2D(hname.c_str(), ("chi2pid vs p: " + std::get<1>(particle_types[i])).c_str(),
                                            nBins, pMin, pMax, nBins, chi2pidMin, chi2pidMax);
            h_mc_chi2pid_vs_p[i]->GetXaxis()->SetTitle("p (GeV)");
            h_mc_chi2pid_vs_p[i]->GetYaxis()->SetTitle("chi2pid");
            h_mc_chi2pid_vs_p[i]->SetStats(false);  // Hide the stat box
        }
    }

    // Create the 1D histograms for each momentum bin
    for (size_t i = 0; i < pBins.size() - 1; ++i) {
        h_data_beta_bins_pos[i] = new TH1D(("h_data_beta_pos_bin_" + std::to_string(i)).c_str(),
                                           ("#beta: " + std::to_string(pBins[i]) + " < p < " + std::to_string(pBins[i + 1])).c_str(),
                                           nBins, 0.7, 1.05);
        h_data_beta_bins_pos[i]->GetXaxis()->SetTitle("#beta");
        h_data_beta_bins_pos[i]->GetYaxis()->SetTitle("Normalized Counts");
        h_data_beta_bins_pos[i]->SetStats(false);

        h_data_beta_bins_neg[i] = new TH1D(("h_data_beta_neg_bin_" + std::to_string(i)).c_str(),
                                           ("#beta: " + std::to_string(pBins[i]) + " < p < " + std::to_string(pBins[i + 1])).c_str(),
                                           nBins, 0.7, 1.05);
        h_data_beta_bins_neg[i]->GetXaxis()->SetTitle("#beta");
        h_data_beta_bins_neg[i]->GetYaxis()->SetTitle("Normalized Counts");
        h_data_beta_bins_neg[i]->SetStats(false);

        if (mcReader) {
            h_mc_beta_bins_pos[i] = new TH1D(("h_mc_beta_pos_bin_" + std::to_string(i)).c_str(),
                                             ("#beta: " + std::to_string(pBins[i]) + " < p < " + std::to_string(pBins[i + 1])).c_str(),
                                             nBins, 0.7, 1.05);
            h_mc_beta_bins_pos[i]->GetXaxis()->SetTitle("#beta");
            h_mc_beta_bins_pos[i]->GetYaxis()->SetTitle("Normalized Counts");
            h_mc_beta_bins_pos[i]->SetLineColor(kRed);
            h_mc_beta_bins_pos[i]->SetStats(false);

            h_mc_beta_bins_neg[i] = new TH1D(("h_mc_beta_neg_bin_" + std::to_string(i)).c_str(),
                                             ("#beta: " + std::to_string(pBins[i]) + " < p < " + std::to_string(pBins[i + 1])).c_str(),
                                             nBins, 0.7, 1.05);
            h_mc_beta_bins_neg[i]->GetXaxis()->SetTitle("#beta");
            h_mc_beta_bins_neg[i]->GetYaxis()->SetTitle("Normalized Counts");
            h_mc_beta_bins_neg[i]->SetLineColor(kRed);
            h_mc_beta_bins_neg[i]->SetStats(false);
        }
    }

    // Fill histograms for data
    while (dataReader.Next()) {
        if (*track_sector_6 != -9999) {  // CD check
            for (size_t i = 0; i < particle_types.size(); ++i) {
                if (*particle_pid == std::get<0>(particle_types[i])) {
                    h_data[i]->Fill(*particle_chi2pid);
                    h_data_chi2pid_vs_p[i]->Fill(*particle_p, *particle_chi2pid);
                    if (*particle_pid > 0) {
                        h_data_beta_vs_p_pos->Fill(*particle_p, *particle_beta);
                        for (size_t bin = 0; bin < pBins.size() - 1; ++bin) {
                            if (*particle_p >= pBins[bin] && *particle_p < pBins[bin + 1]) {
                                h_data_beta_bins_pos[bin]->Fill(*particle_beta);
                                break;
                            }
                        }
                    } else {
                        h_data_beta_vs_p_neg->Fill(*particle_p, *particle_beta);
                        for (size_t bin = 0; bin < pBins.size() - 1; ++bin) {
                            if (*particle_p >= pBins[bin] && *particle_p < pBins[bin + 1]) {
                                h_data_beta_bins_neg[bin]->Fill(*particle_beta);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    // Fill histograms for MC (if applicable)
    if (mcReader) {
        while (mcReader->Next()) {
            if (**mc_track_sector_6 != -9999) {  // CD check
                for (size_t i = 0; i < particle_types.size(); ++i) {
                    if (**mc_particle_pid == std::get<0>(particle_types[i])) {
                        h_mc[i]->Fill(**mc_particle_chi2pid);
                        h_mc_chi2pid_vs_p[i]->Fill(**mc_particle_p, **mc_particle_chi2pid);
                        if (**mc_particle_pid > 0) {
                            h_mc_beta_vs_p_pos->Fill(**mc_particle_p, **mc_particle_beta);
                            for (size_t bin = 0; bin < pBins.size() - 1; ++bin) {
                                if (**mc_particle_p >= pBins[bin] && **mc_particle_p < pBins[bin + 1]) {
                                    h_mc_beta_bins_pos[bin]->Fill(**mc_particle_beta);
                                    break;
                                }
                            }
                        } else {
                            h_mc_beta_vs_p_neg->Fill(**mc_particle_p, **mc_particle_beta);
                            for (size_t bin = 0; bin < pBins.size() - 1; ++bin) {
                                if (**mc_particle_p >= pBins[bin] && **mc_particle_p < pBins[bin + 1]) {
                                    h_mc_beta_bins_neg[bin]->Fill(**mc_particle_beta);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Normalize 1D histograms and set y-axis range
    double max_y = 0;
    for (size_t i = 0; i < particle_types.size(); ++i) {
        h_data[i]->Scale(1.0 / h_data[i]->Integral());
        if (h_data[i]->GetMaximum() > max_y) max_y = h_data[i]->GetMaximum();

        if (mcReader) {
            h_mc[i]->Scale(1.0 / h_mc[i]->Integral());
            if (h_mc[i]->GetMaximum() > max_y) max_y = h_mc[i]->GetMaximum();
        }
    }

    // Update y-axis range to be 1.2 times the maximum value
    for (size_t i = 0; i < particle_types.size(); ++i) {
        h_data[i]->SetMaximum(1.2 * max_y);
        if (mcReader) h_mc[i]->SetMaximum(1.2 * max_y);
    }

    // Normalize beta histograms for momentum bins and find max_y for each bin
    for (size_t bin = 0; bin < pBins.size() - 1; ++bin) {
        double max_y_bin = 0;

        h_data_beta_bins_pos[bin]->Scale(1.0 / h_data_beta_bins_pos[bin]->Integral());
        h_data_beta_bins_neg[bin]->Scale(1.0 / h_data_beta_bins_neg[bin]->Integral());

        if (h_data_beta_bins_pos[bin]->GetMaximum() > max_y_bin) max_y_bin = h_data_beta_bins_pos[bin]->GetMaximum();
        if (h_data_beta_bins_neg[bin]->GetMaximum() > max_y_bin) max_y_bin = h_data_beta_bins_neg[bin]->GetMaximum();

        if (mcReader) {
            h_mc_beta_bins_pos[bin]->Scale(1.0 / h_mc_beta_bins_pos[bin]->Integral());
            h_mc_beta_bins_neg[bin]->Scale(1.0 / h_mc_beta_bins_neg[bin]->Integral());

            if (h_mc_beta_bins_pos[bin]->GetMaximum() > max_y_bin) max_y_bin = h_mc_beta_bins_pos[bin]->GetMaximum();
            if (h_mc_beta_bins_neg[bin]->GetMaximum() > max_y_bin) max_y_bin = h_mc_beta_bins_neg[bin]->GetMaximum();
        }

        h_data_beta_bins_pos[bin]->SetMaximum(1.2 * max_y_bin);
        h_data_beta_bins_neg[bin]->SetMaximum(1.2 * max_y_bin);

        if (mcReader) {
            h_mc_beta_bins_pos[bin]->SetMaximum(1.2 * max_y_bin);
            h_mc_beta_bins_neg[bin]->SetMaximum(1.2 * max_y_bin);
        }
    }

    // Draw 1D histograms on the same canvas for each particle type
    for (size_t i = 0; i < particle_types.size(); ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.15);  // Add padding to the left
        h_data[i]->Draw("HIST");
        if (mcReader) h_mc[i]->Draw("HIST SAME");

            TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->AddEntry(h_data[i], "Data", "l");
            if (mcReader) legend->AddEntry(h_mc[i], "MC", "l");
            // legend->AddEntry(h_data[i], "NH_{3}", "l");
            // if (mcReader) legend->AddEntry(h_mc[i], "C", "l");
            // legend->AddEntry(h_data[i], "NH_{3} pass-1", "l");
            // if (mcReader) legend->AddEntry(h_mc[i], "NH_{3} preliminary", "l");
            legend->Draw();
        }

        // Draw 2D histograms for data
        for (size_t i = 0; i < particle_types.size(); ++i) {
            c_data_2D->cd(i + 1);
            gPad->SetLeftMargin(0.15);  // Add padding to the left
            gPad->SetLogz();  // Set log scale for z-axis
            h_data_chi2pid_vs_p[i]->Draw("COLZ");
        }

        // Draw 2D histograms for MC (if applicable)
        if (mcReader) {
            for (size_t i = 0; i < particle_types.size(); ++i) {
                c_mc_2D->cd(i + 1);
                gPad->SetLeftMargin(0.15);  // Add padding to the left
                gPad->SetLogz();  // Set log scale for z-axis
                h_mc_chi2pid_vs_p[i]->Draw("COLZ");
            }
        }

        // Draw the combined beta vs p histograms
        c_data_pos_neg_beta->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetLogz();
        h_data_beta_vs_p_pos->Draw("COLZ");

        c_data_pos_neg_beta->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetLogz();
        h_data_beta_vs_p_neg->Draw("COLZ");

        if (mcReader) {
            c_mc_pos_neg_beta->cd(1);
            gPad->SetLeftMargin(0.15);
            gPad->SetLogz();
            h_mc_beta_vs_p_pos->Draw("COLZ");

            c_mc_pos_neg_beta->cd(2);
            gPad->SetLeftMargin(0.15);
            gPad->SetLogz();
            h_mc_beta_vs_p_neg->Draw("COLZ");
        }

        // Save the canvases
        c->SaveAs("output/calibration/cvt/chi2pid/chi2pid_cd.png");
        c_data_2D->SaveAs("output/calibration/cvt/chi2pid/chi2pid_vs_p_cd.png");
        c_data_pos_neg_beta->SaveAs("output/calibration/cvt/chi2pid/beta_vs_p_cd.png");

        if (mcReader) {
            c_mc_2D->SaveAs("output/calibration/cvt/chi2pid/mc_chi2pid_vs_p_cd.png");
            c_mc_pos_neg_beta->SaveAs("output/calibration/cvt/chi2pid/mc_beta_vs_p_cd.png");
        }

        // Create 4x4 canvases for beta vs p momentum bin histograms
        TCanvas* c_data_beta_bins_pos = new TCanvas("c_data_beta_bins_pos", "Beta vs p Binned (Positive Tracks)", 1800, 1800);
        c_data_beta_bins_pos->Divide(4, 4);
        TCanvas* c_data_beta_bins_neg = new TCanvas("c_data_beta_bins_neg", "Beta vs p Binned (Negative Tracks)", 1800, 1800);
        c_data_beta_bins_neg->Divide(4, 4);

        TCanvas* c_mc_beta_bins_pos = nullptr;
        TCanvas* c_mc_beta_bins_neg = nullptr;

        if (mcReader) {
            c_mc_beta_bins_pos = new TCanvas("c_mc_beta_bins_pos", "MC Beta vs p Binned (Positive Tracks)", 1800, 1800);
            c_mc_beta_bins_pos->Divide(4, 4);
            c_mc_beta_bins_neg = new TCanvas("c_mc_beta_bins_neg", "MC Beta vs p Binned (Negative Tracks)", 1800, 1800);
            c_mc_beta_bins_neg->Divide(4, 4);
        }

        // Fill 4x4 canvases for beta bins
        for (size_t bin = 0; bin < pBins.size() - 1; ++bin) {
            c_data_beta_bins_pos->cd(bin + 1);
            gPad->SetLeftMargin(0.15);
            h_data_beta_bins_pos[bin]->Draw("HIST");
            if (mcReader) {
                h_mc_beta_bins_pos[bin]->Draw("HIST SAME");
            }
            TLegend* legend_pos = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend_pos->AddEntry(h_data_beta_bins_pos[bin], "Data", "l");
            if (mcReader) legend_pos->AddEntry(h_mc_beta_bins_pos[bin], "MC", "l");
            // legend_pos->AddEntry(h_data_beta_bins_pos[bin], "NH_{3} pass-1", "l");
            // if (mcReader) legend_pos->AddEntry(h_mc_beta_bins_pos[bin], "C pass-1", "l");
            legend_pos->Draw();

            c_data_beta_bins_neg->cd(bin + 1);
            gPad->SetLeftMargin(0.15);
            h_data_beta_bins_neg[bin]->Draw("HIST");
            if (mcReader) {
                h_mc_beta_bins_neg[bin]->Draw("HIST SAME");
            }
            TLegend* legend_neg = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend_neg->AddEntry(h_data_beta_bins_neg[bin], "Data", "l");
            if (mcReader) legend_neg->AddEntry(h_mc_beta_bins_neg[bin], "MC", "l");
            legend_neg->Draw();
        }

        // Save the 4x4 canvases
        c_data_beta_bins_pos->SaveAs("output/calibration/cvt/chi2pid/beta_vs_p_binned_pos_cd.png");
        c_data_beta_bins_neg->SaveAs("output/calibration/cvt/chi2pid/beta_vs_p_binned_neg_cd.png");

        // Clean up
        delete c;
        delete c_data_2D;
        delete c_data_pos_neg_beta;
        delete c_data_beta_bins_pos;
        delete c_data_beta_bins_neg;

        if (mcReader) {
            delete c_mc_2D;
            delete c_mc_pos_neg_beta;
            delete c_mc_beta_bins_pos;
            delete c_mc_beta_bins_neg;
        }

        for (size_t i = 0; i < particle_types.size(); ++i) {
            delete h_data[i];
            delete h_data_chi2pid_vs_p[i];
            if (mcReader) {
                delete h_mc[i];
                delete h_mc_chi2pid_vs_p[i];
            }
        }

        for (size_t bin = 0; bin < pBins.size() - 1; ++bin) {
            delete h_data_beta_bins_pos[bin];
            delete h_data_beta_bins_neg[bin];
            if (mcReader) {
                delete h_mc_beta_bins_pos[bin];
                delete h_mc_beta_bins_neg[bin];
            }
        }

        if (mc_particle_chi2pid) delete mc_particle_chi2pid;
        if (mc_particle_p) delete mc_particle_p;
        if (mc_particle_beta) delete mc_particle_beta;
        if (mc_track_sector_6) delete mc_track_sector_6;
        if (mc_particle_pid) delete mc_particle_pid;
    }


// Helper functions for theta_CVT and phi_CVT calculations
double calculate_phi(double x, double y) {
    double phi = atan2(x, y) * 180.0 / M_PI;
    phi = phi - 90;
    if (phi < 0) {
        phi += 360;
    }
    phi = 360 - phi;
    return phi;
}

double calculate_theta(double x, double y, double z) {
    double r = sqrt(x*x + y*y + z*z);
    return acos(z / r) * 180.0 / M_PI;
}

bool cvt_fiducial(double edge_1, double edge_3, double edge_5, double edge_7, 
     double edge_12) {
    // Dummy function, always returns true
    return edge_1 > 0 && edge_3 > 0 && edge_5 > 0 && edge_7 > -2 && edge_12 > -5;
    return true;
}

void plot_chi2_ndf_vs_phi_CVT_2D(TTreeReader& dataReader, TTreeReader* mcReader, const std::vector<std::tuple<int, std::string, std::string>>& particle_types) {
    int nBins = 100;

    // Declare TTreeReaderValues for trajectory and track variables for data
    TTreeReaderValue<double> traj_x_12(dataReader, "traj_x_12");
    TTreeReaderValue<double> traj_y_12(dataReader, "traj_y_12");
    TTreeReaderValue<double> traj_z_12(dataReader, "traj_z_12");
    TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");
    TTreeReaderValue<double> track_chi2_5(dataReader, "track_chi2_5");
    TTreeReaderValue<int> track_ndf_5(dataReader, "track_ndf_5");

    // Declare TTreeReaderValues for trajectory and track variables for MC if available
    TTreeReaderValue<double>* mc_traj_x_12 = nullptr;
    TTreeReaderValue<double>* mc_traj_y_12 = nullptr;
    TTreeReaderValue<double>* mc_traj_z_12 = nullptr;
    TTreeReaderValue<int>* mc_particle_pid = nullptr;
    TTreeReaderValue<double>* mc_track_chi2_5 = nullptr;
    TTreeReaderValue<int>* mc_track_ndf_5 = nullptr;

    if (mcReader) {
        mc_traj_x_12 = new TTreeReaderValue<double>(*mcReader, "traj_x_12");
        mc_traj_y_12 = new TTreeReaderValue<double>(*mcReader, "traj_y_12");
        mc_traj_z_12 = new TTreeReaderValue<double>(*mcReader, "traj_z_12");
        mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
        mc_track_chi2_5 = new TTreeReaderValue<double>(*mcReader, "track_chi2_5");
        mc_track_ndf_5 = new TTreeReaderValue<int>(*mcReader, "track_ndf_5");
    }

    for (const auto& particle_type : particle_types) {
        int pid = std::get<0>(particle_type);
        std::string particle_name = std::get<1>(particle_type);
        std::string particle_latex = std::get<2>(particle_type);

        // Create histograms for chi2/ndf vs phi_CVT for data and MC
        TH2D* h_chi2_vs_phi_CVT_data = new TH2D(("h_chi2_vs_phi_CVT_data_" + particle_name).c_str(), 
                                                ("#chi^{2}/ndf vs #phi_{CVT} (Data, " + particle_latex + ")").c_str(), 
                                                nBins, 0, 360, nBins, 0, 100);
        h_chi2_vs_phi_CVT_data->GetXaxis()->SetTitle("#phi_{CVT}");
        h_chi2_vs_phi_CVT_data->GetYaxis()->SetTitle("#chi^{2}/ndf");

        TH2D* h_chi2_vs_phi_CVT_mc = nullptr;
        if (mcReader) {
            h_chi2_vs_phi_CVT_mc = new TH2D(("h_chi2_vs_phi_CVT_mc_" + particle_name).c_str(), 
                                            ("#chi^{2}/ndf vs #phi_{CVT} (MC, " + particle_latex + ")").c_str(), 
                                            nBins, 0, 360, nBins, 0, 100);
            h_chi2_vs_phi_CVT_mc->GetXaxis()->SetTitle("#phi_{CVT}");
            h_chi2_vs_phi_CVT_mc->GetYaxis()->SetTitle("#chi^{2}/ndf");
        }

        // Fill the histograms for data
        dataReader.Restart();
        while (dataReader.Next()) {
            if (*particle_pid == pid && *track_ndf_5 > 0 && *track_chi2_5 < 10000) {
                double chi2_ndf = *track_chi2_5 / *track_ndf_5;

                if (*traj_x_12 != -9999 && *traj_y_12 != -9999 && *traj_z_12 != -9999) {
                    double phi_CVT = calculate_phi(*traj_x_12, *traj_y_12);
                    // double theta_CVT = calculate_theta(*traj_x_12, *traj_y_12, *traj_z_12);

                    h_chi2_vs_phi_CVT_data->Fill(phi_CVT, chi2_ndf);
                }
            }
        }

        // Fill the histograms for MC if available
        if (mcReader) {
            mcReader->Restart();
            while (mcReader->Next()) {
                if (**mc_particle_pid == pid && **mc_track_ndf_5 > 0 && **mc_track_chi2_5 < 10000) {
                    double mc_chi2_ndf = **mc_track_chi2_5 / **mc_track_ndf_5;

                    if (**mc_traj_x_12 != -9999 && **mc_traj_y_12 != -9999 && **mc_traj_z_12 != -9999) {
                        double mc_phi_CVT = calculate_phi(**mc_traj_x_12, **mc_traj_y_12);
                        // double mc_theta_CVT = calculate_theta(**mc_traj_x_12, **mc_traj_y_12, **mc_traj_z_12);

                        h_chi2_vs_phi_CVT_mc->Fill(mc_phi_CVT, mc_chi2_ndf);
                    }
                }
            }
        }

        // Save the histograms
        TCanvas* c_data = new TCanvas(("c_chi2_vs_phi_CVT_data_" + particle_name).c_str(), ("#chi^{2}/ndf vs #phi_{CVT} (Data, " + particle_latex + ")").c_str(), 800, 600);
        gPad->SetLogz();  // Set the y-axis to log scale
        h_chi2_vs_phi_CVT_data->Draw("COLZ");
        c_data->SaveAs(("output/calibration/cvt/determination/chi2_vs_phi_CVT_data_" + particle_name + ".png").c_str());
        delete c_data;

        if (mcReader) {
            TCanvas* c_mc = new TCanvas(("c_chi2_vs_phi_CVT_mc_" + particle_name).c_str(), ("#chi^{2}/ndf vs #phi_{CVT} (MC, " + particle_latex + ")").c_str(), 800, 600);
            gPad->SetLogz();  // Set the y-axis to log scale
            h_chi2_vs_phi_CVT_mc->Draw("COLZ");
            c_mc->SaveAs(("output/calibration/cvt/determination/chi2_vs_phi_CVT_mc_" + particle_name + ".png").c_str());
            delete c_mc;
            delete h_chi2_vs_phi_CVT_mc;
        }

        delete h_chi2_vs_phi_CVT_data;
    }

    // Clean up dynamically allocated memory for MC
    if (mcReader) {
        delete mc_traj_x_12;
        delete mc_traj_y_12;
        delete mc_traj_z_12;
        delete mc_particle_pid;
        delete mc_track_chi2_5;
        delete mc_track_ndf_5;
    }
}

void cvt_fiducial_determination(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    int nBinsX = 50;  // Number of bins for edge
    int nBinsY = 100;  // Number of bins for chi2/ndf

    // Define CVT layers and edges for data
    std::vector<std::tuple<TTreeReaderValue<double>*, std::string, double, double>> layers = {
        {new TTreeReaderValue<double>(dataReader, "traj_edge_1"), "layer_1", -2, 2.2},
        {new TTreeReaderValue<double>(dataReader, "traj_edge_3"), "layer_3", -2, 2.2},
        {new TTreeReaderValue<double>(dataReader, "traj_edge_5"), "layer_5", -2, 2.2},
        {new TTreeReaderValue<double>(dataReader, "traj_edge_7"), "layer_7", -5, 15},
        {new TTreeReaderValue<double>(dataReader, "traj_edge_12"), "layer_12", -10, 25}
    };

    // Define CVT layers and edges for MC if available
    std::vector<std::tuple<TTreeReaderValue<double>*, std::string, double, double>> mc_layers;
    if (mcReader) {
        mc_layers = {
            {new TTreeReaderValue<double>(*mcReader, "traj_edge_1"), "layer_1", -2, 2.2},
            {new TTreeReaderValue<double>(*mcReader, "traj_edge_3"), "layer_3", -2, 2.2},
            {new TTreeReaderValue<double>(*mcReader, "traj_edge_5"), "layer_5", -2, 2.2},
            {new TTreeReaderValue<double>(*mcReader, "traj_edge_7"), "layer_7", -5, 15},
            {new TTreeReaderValue<double>(*mcReader, "traj_edge_12"), "layer_12", -10, 25}
        };
    }

    // Define particle PID and track variables for data
    TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");
    TTreeReaderValue<double> track_chi2_5(dataReader, "track_chi2_5");
    TTreeReaderValue<int> track_ndf_5(dataReader, "track_ndf_5");

    // Define particle PID and track variables for MC if available
    TTreeReaderValue<int>* mc_particle_pid = nullptr;
    TTreeReaderValue<double>* mc_track_chi2_5 = nullptr;
    TTreeReaderValue<int>* mc_track_ndf_5 = nullptr;

    if (mcReader) {
        mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
        mc_track_chi2_5 = new TTreeReaderValue<double>(*mcReader, "track_chi2_5");
        mc_track_ndf_5 = new TTreeReaderValue<int>(*mcReader, "track_ndf_5");
    }

    std::vector<std::tuple<int, std::string, std::string>> particle_types = {
        {211, "pip", "#pi^{+}"},
        {-211, "pim", "#pi^{-}"},
        {321, "kp", "k^{+}"},
        {-321, "km", "k^{-}"},
        {2212, "proton", "proton"}
    };

    plot_chi2_ndf_vs_phi_CVT_2D(dataReader, mcReader, particle_types);

    for (const auto& particle_type : particle_types) {
        int pid = std::get<0>(particle_type);
        std::string particle_name = std::get<1>(particle_type);
        std::string particle_latex = std::get<2>(particle_type);

        std::vector<TH2D*> h_chi2_vs_edge_data, h_chi2_vs_edge_mc;

        // Initialize 2D histograms for data and MC
        for (const auto& layer : layers) {
            std::string layer_name = std::get<1>(layer);
            double xMin = std::get<2>(layer);
            double xMax = std::get<3>(layer);
            double yMin = 0;
            double yMax = 100;  // This can be adjusted depending on your expected range of chi2/ndf

            h_chi2_vs_edge_data.push_back(new TH2D(("h_chi2_vs_edge_data_" + layer_name + "_" + particle_name).c_str(), (particle_latex + " - " + layer_name).c_str(), nBinsX, xMin, xMax, nBinsY, yMin, yMax));

            if (mcReader) {
                h_chi2_vs_edge_mc.push_back(new TH2D(("h_chi2_vs_edge_mc_" + layer_name + "_" + particle_name).c_str(), (particle_latex + " - MC - " + layer_name).c_str(), nBinsX, xMin, xMax, nBinsY, yMin, yMax));
            }
        }

        // Fill 2D histograms for data
        dataReader.Restart();
        while (dataReader.Next()) {
            if (*particle_pid == pid && *track_ndf_5 > 0 && *track_chi2_5 < 10000) {
                double chi2_ndf = *track_chi2_5 / *track_ndf_5;

                for (size_t i = 0; i < layers.size(); ++i) {
                    double traj_edge = **std::get<0>(layers[i]);

                    if (traj_edge != -9999) {
                        h_chi2_vs_edge_data[i]->Fill(traj_edge, chi2_ndf);
                    }
                }
            }
        }

        // Fill 2D histograms for MC
        if (mcReader) {
            mcReader->Restart();
            while (mcReader->Next()) {
                if (**mc_particle_pid == pid && **mc_track_ndf_5 > 0 && **mc_track_chi2_5 < 10000) {
                    double mc_chi2_ndf = **mc_track_chi2_5 / **mc_track_ndf_5;

                    for (size_t i = 0; i < mc_layers.size(); ++i) {
                        double traj_edge = **std::get<0>(mc_layers[i]);

                        if (traj_edge != -9999) {
                            h_chi2_vs_edge_mc[i]->Fill(traj_edge, mc_chi2_ndf);
                        }
                    }
                }
            }
        }

        // Plot data results
        TCanvas* c_layer_data = new TCanvas(("c_data_" + particle_name).c_str(), ("#chi^{2}/ndf vs edge for " + particle_latex + " (Data)").c_str(), 1800, 1200);
        c_layer_data->Divide(3, 2);
        c_layer_data->SetLeftMargin(0.2);  // Increase left margin for better y-axis label visibility

        for (size_t i = 0; i < layers.size(); ++i) {
            c_layer_data->cd(i + 1);
            gPad->SetLogz();  // Set the z-axis to log scale
            gPad->SetMargin(0.2, 0.15, 0.2, 0.1);  // Adjust margins: left, right, bottom, top
            h_chi2_vs_edge_data[i]->GetXaxis()->SetTitle("edge (cm)");
            h_chi2_vs_edge_data[i]->GetYaxis()->SetTitle("#chi^{2}/ndf");
            h_chi2_vs_edge_data[i]->Draw("COLZ");
        }

        c_layer_data->SaveAs(("output/calibration/cvt/determination/chi2_per_ndf_vs_edge_data_" + particle_name + ".png").c_str());

        // Plot MC results
        if (mcReader) {
            TCanvas* c_layer_mc = new TCanvas(("c_mc_" + particle_name).c_str(), ("#chi^{2}/ndf vs edge for " + particle_latex + " (MC)").c_str(), 1800, 1200);
            c_layer_mc->Divide(3, 2);
            c_layer_mc->SetLeftMargin(0.2);  // Increase left margin for better y-axis label visibility

            for (size_t i = 0; i < mc_layers.size(); ++i) {
                c_layer_mc->cd(i + 1);
                gPad->SetLogz();  // Set the z-axis to log scale
                gPad->SetMargin(0.2, 0.15, 0.2, 0.1);  // Adjust margins: left, right, bottom, top
                h_chi2_vs_edge_mc[i]->GetXaxis()->SetTitle("edge (cm)");
                h_chi2_vs_edge_mc[i]->GetYaxis()->SetTitle("#chi^{2}/ndf");
                h_chi2_vs_edge_mc[i]->Draw("COLZ");
            }

            c_layer_mc->SaveAs(("output/calibration/cvt/determination/chi2_per_ndf_vs_edge_mc_" + particle_name + ".png").c_str());
            delete c_layer_mc;
        }
        // Clean up
        for (auto hist : h_chi2_vs_edge_data) delete hist;
        if (mcReader) {
            for (auto hist : h_chi2_vs_edge_mc) delete hist;
        }
        delete c_layer_data;
    }

    // Clean up dynamically allocated memory for MC
    if (mcReader) {
        for (auto& mc_layer : mc_layers) {
            delete std::get<0>(mc_layer);
        }

        delete mc_particle_pid;
        delete mc_track_chi2_5;
        delete mc_track_ndf_5;
    }

    // Clean up dynamically allocated memory for data
    for (auto& layer : layers) {
        delete std::get<0>(layer);
    }
}

void plot_cvt_hit_position(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
    int nBins = 100;

    std::vector<std::tuple<std::string, std::string, std::string, double, double>> layers = {
        {"traj_x_1", "traj_y_1", "layer_1", -10, 10},
        {"traj_x_3", "traj_y_3", "layer_3", -12, 12},
        {"traj_x_5", "traj_y_5", "layer_5", -15, 15},
        {"traj_x_7", "traj_y_7", "layer_7", -17.5, 17.5},
        {"traj_x_12", "traj_y_12", "layer_12", -25, 25}
    };

    std::vector<std::tuple<int, std::string, std::string>> particle_types = {
        {211, "pip", "#pi^{+}"},
        {-211, "pim", "#pi^{-}"},
        {321, "kp", "k^{+}"},
        {-321, "km", "k^{-}"},
        {2212, "proton", "proton"}
    };

    // Declare TTreeReaderValues for the CVT edge and track variables
    TTreeReaderValue<double> traj_edge_1(dataReader, "traj_edge_1");
    TTreeReaderValue<double> traj_edge_3(dataReader, "traj_edge_3");
    TTreeReaderValue<double> traj_edge_5(dataReader, "traj_edge_5");
    TTreeReaderValue<double> traj_edge_7(dataReader, "traj_edge_7");
    TTreeReaderValue<double> traj_edge_12(dataReader, "traj_edge_12");

    TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");
    TTreeReaderValue<double> theta(dataReader, "theta");

    TTreeReaderValue<double> traj_x_12(dataReader, "traj_x_12");
    TTreeReaderValue<double> traj_y_12(dataReader, "traj_y_12");
    TTreeReaderValue<double> traj_z_12(dataReader, "traj_z_12");

    TTreeReaderValue<double>* mc_traj_edge_1 = nullptr;
    TTreeReaderValue<double>* mc_traj_edge_3 = nullptr;
    TTreeReaderValue<double>* mc_traj_edge_5 = nullptr;
    TTreeReaderValue<double>* mc_traj_edge_7 = nullptr;
    TTreeReaderValue<double>* mc_traj_edge_12 = nullptr;
    TTreeReaderValue<int>* mc_particle_pid = nullptr;
    TTreeReaderValue<double>* mc_theta = nullptr;
    TTreeReaderValue<double>* mc_traj_x_12 = nullptr;
    TTreeReaderValue<double>* mc_traj_y_12 = nullptr;
    TTreeReaderValue<double>* mc_traj_z_12 = nullptr;

    if (mcReader) {
        mc_traj_edge_1 = new TTreeReaderValue<double>(*mcReader, "traj_edge_1");
        mc_traj_edge_3 = new TTreeReaderValue<double>(*mcReader, "traj_edge_3");
        mc_traj_edge_5 = new TTreeReaderValue<double>(*mcReader, "traj_edge_5");
        mc_traj_edge_7 = new TTreeReaderValue<double>(*mcReader, "traj_edge_7");
        mc_traj_edge_12 = new TTreeReaderValue<double>(*mcReader, "traj_edge_12");
        mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
        mc_theta = new TTreeReaderValue<double>(*mcReader, "theta");
        mc_traj_x_12 = new TTreeReaderValue<double>(*mcReader, "traj_x_12");
        mc_traj_y_12 = new TTreeReaderValue<double>(*mcReader, "traj_y_12");
        mc_traj_z_12 = new TTreeReaderValue<double>(*mcReader, "traj_z_12");
    }

    // Declare TTreeReaderValues for trajectory x and y coordinates
    std::vector<TTreeReaderValue<double>> traj_x;
    std::vector<TTreeReaderValue<double>> traj_y;

    std::vector<TTreeReaderValue<double>*> mc_traj_x;
    std::vector<TTreeReaderValue<double>*> mc_traj_y;

    // Initialize TTreeReaderValues for each layer
    for (const auto& layer : layers) {
        traj_x.emplace_back(dataReader, std::get<0>(layer).c_str());
        traj_y.emplace_back(dataReader, std::get<1>(layer).c_str());

        if (mcReader) {
            mc_traj_x.push_back(new TTreeReaderValue<double>(*mcReader, std::get<0>(layer).c_str()));
            mc_traj_y.push_back(new TTreeReaderValue<double>(*mcReader, std::get<1>(layer).c_str()));
        }
    }

    for (const auto& particle_type : particle_types) {
        int pid = std::get<0>(particle_type);
        std::string particle_name = std::get<1>(particle_type);
        std::string particle_latex = std::get<2>(particle_type);

        // Create a canvas for data
        TCanvas* c_data = new TCanvas(("c_data_" + particle_name).c_str(), ("Data CVT Hit Position (" + particle_latex + ")").c_str(), 1800, 1200);
        c_data->Divide(5, 2);

        TCanvas* c_mc = nullptr;
        if (mcReader) {
            c_mc = new TCanvas(("c_mc_" + particle_name).c_str(), ("MC CVT Hit Position (" + particle_latex + ")").c_str(), 1800, 1200);
            c_mc->Divide(5, 2);
        }

        // Create histograms for data and MC
        std::vector<TH2D*> h_data_before(5), h_data_after(5);
        std::vector<TH2D*> h_mc_before(5), h_mc_after(5);

        for (int layer_idx = 0; layer_idx < 5; ++layer_idx) {
            std::string layer_name = std::get<2>(layers[layer_idx]);
            double xMin = std::get<3>(layers[layer_idx]);
            double xMax = std::get<4>(layers[layer_idx]);
            double yMin = xMin;
            double yMax = xMax;

            h_data_before[layer_idx] = new TH2D(("h_data_before_" + layer_name).c_str(), ("Data " + layer_name + " Before Cuts (" + particle_latex + ")").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
            h_data_after[layer_idx] = new TH2D(("h_data_after_" + layer_name).c_str(), ("Data " + layer_name + " After Cuts (" + particle_latex + ")").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
            h_data_before[layer_idx]->GetXaxis()->SetTitle("x");
            h_data_before[layer_idx]->GetYaxis()->SetTitle("y");
            h_data_after[layer_idx]->GetXaxis()->SetTitle("x");
            h_data_after[layer_idx]->GetYaxis()->SetTitle("y");

            if (mcReader) {
                h_mc_before[layer_idx] = new TH2D(("h_mc_before_" + layer_name).c_str(), ("MC " + layer_name + " Before Cuts (" + particle_latex + ")").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
                h_mc_after[layer_idx] = new TH2D(("h_mc_after_" + layer_name).c_str(), ("MC " + layer_name + " After Cuts (" + particle_latex + ")").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
                h_mc_before[layer_idx]->GetXaxis()->SetTitle("x");
                h_mc_before[layer_idx]->GetYaxis()->SetTitle("y");
                h_mc_after[layer_idx]->GetXaxis()->SetTitle("x");
                h_mc_after[layer_idx]->GetYaxis()->SetTitle("y");
            }
        }

        // Fill the data histograms
        dataReader.Restart();
        while (dataReader.Next()) {
            if (*particle_pid == pid) {
                for (int layer_idx = 0; layer_idx < 5; ++layer_idx) {
                    double traj_x_value = *traj_x[layer_idx];
                    double traj_y_value = *traj_y[layer_idx];

                    if (traj_x_value != -9999 && traj_y_value != -9999) {
                        h_data_before[layer_idx]->Fill(traj_x_value, traj_y_value);
                        if (cvt_fiducial(*traj_edge_1, *traj_edge_3, *traj_edge_5, *traj_edge_7, *traj_edge_12)) {
                            h_data_after[layer_idx]->Fill(traj_x_value, traj_y_value);
                        }
                    }
                }
            }
        }

        // Fill the MC histograms if available
        if (mcReader) {
            mcReader->Restart();
            while (mcReader->Next()) {
                if (**mc_particle_pid == pid) {
                    for (int layer_idx = 0; layer_idx < 5; ++layer_idx) {
                        double mc_traj_x_value = **mc_traj_x[layer_idx];
                        double mc_traj_y_value = **mc_traj_y[layer_idx];
                        if (mc_traj_x_value != -9999 && mc_traj_y_value != -9999) {
                            h_mc_before[layer_idx]->Fill(mc_traj_x_value, mc_traj_y_value);
                            if (cvt_fiducial(**mc_traj_edge_1, **mc_traj_edge_3, **mc_traj_edge_5, **mc_traj_edge_7, **mc_traj_edge_12)) {
                                h_mc_after[layer_idx]->Fill(mc_traj_x_value, mc_traj_y_value);
                            }
                        }
                    }
                }
            }
        }
        // Find the maximum value across all histograms for consistent scaling
        double max_value_data = 0, max_value_mc = 0;
        for (int layer_idx = 0; layer_idx < 5; ++layer_idx) {
            max_value_data = std::max(max_value_data, h_data_before[layer_idx]->GetMaximum());
            max_value_data = std::max(max_value_data, h_data_after[layer_idx]->GetMaximum());
            if (mcReader) {
                max_value_mc = std::max(max_value_mc, h_mc_before[layer_idx]->GetMaximum());
                max_value_mc = std::max(max_value_mc, h_mc_after[layer_idx]->GetMaximum());
            }
        }

        // Set the maximum for each histogram to ensure consistent scaling
        for (int layer_idx = 0; layer_idx < 5; ++layer_idx) {
            h_data_before[layer_idx]->SetMaximum(max_value_data * 1.1);
            h_data_after[layer_idx]->SetMaximum(max_value_data * 1.1);
            if (mcReader) {
                h_mc_before[layer_idx]->SetMaximum(max_value_mc * 1.1);
                h_mc_after[layer_idx]->SetMaximum(max_value_mc * 1.1);
            }
        }

        // Draw and save the data canvas
        for (int layer_idx = 0; layer_idx < 5; ++layer_idx) {
            c_data->cd(layer_idx + 1);
            gPad->SetLogz();
            gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
            h_data_before[layer_idx]->Draw("COLZ");

            c_data->cd(layer_idx + 6);
            gPad->SetLogz();
            gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
            h_data_after[layer_idx]->Draw("COLZ");
        }
        c_data->SaveAs(("output/calibration/cvt/positions/data_" + particle_name + "_cvt_hit_position.png").c_str());

        // Draw and save the MC canvas if available
        if (mcReader) {
            for (int layer_idx = 0; layer_idx < 5; ++layer_idx) {
                c_mc->cd(layer_idx + 1);
                gPad->SetLogz();
                gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
                h_mc_before[layer_idx]->Draw("COLZ");

                c_mc->cd(layer_idx + 6);
                gPad->SetLogz();
                gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
                h_mc_after[layer_idx]->Draw("COLZ");
            }
            c_mc->SaveAs(("output/calibration/cvt/positions/mc_" + particle_name + "_cvt_hit_position.png").c_str());
        }

        // Convert theta to degrees for both data and MC
        std::vector<double> theta_in_degrees_data, theta_in_degrees_mc;
        std::vector<double> theta_CVT_data, phi_CVT_data;
        std::vector<double> traj_edge_1_data, traj_edge_3_data, traj_edge_5_data, traj_edge_7_data, traj_edge_12_data;
        std::vector<double> theta_CVT_mc, phi_CVT_mc;
        std::vector<double> traj_edge_1_mc, traj_edge_3_mc, traj_edge_5_mc, traj_edge_7_mc, traj_edge_12_mc;

        // Calculate theta_CVT and phi_CVT for data
        dataReader.Restart();
        while (dataReader.Next()) {
            if (*particle_pid == pid) {
                if (*traj_x_12 != -9999 && *traj_y_12 != -9999 && *traj_z_12 != -9999) {
                    double theta_CVT_value = calculate_theta(*traj_x_12, *traj_y_12, *traj_z_12);
                    double phi_CVT_value = calculate_phi(*traj_x_12, *traj_y_12);

                    // Store the CVT angles and corresponding edge values
                    theta_CVT_data.push_back(theta_CVT_value);
                    phi_CVT_data.push_back(phi_CVT_value);
                    theta_in_degrees_data.push_back(*theta);

                    // Store the corresponding edge values
                    traj_edge_1_data.push_back(*traj_edge_1);
                    traj_edge_3_data.push_back(*traj_edge_3);
                    traj_edge_5_data.push_back(*traj_edge_5);
                    traj_edge_7_data.push_back(*traj_edge_7);
                    traj_edge_12_data.push_back(*traj_edge_12);
                }
            }
        }

        // Calculate theta_CVT and phi_CVT for MC
        if (mcReader) {
            mcReader->Restart();
            while (mcReader->Next()) {
                if (**mc_particle_pid == pid) {
                    if (**mc_traj_x_12 != -9999 && **mc_traj_y_12 != -9999 && **mc_traj_z_12 != -9999) {
                        double mc_theta_CVT_value = calculate_theta(**mc_traj_x_12, **mc_traj_y_12, **mc_traj_z_12);
                        double mc_phi_CVT_value = calculate_phi(**mc_traj_x_12, **mc_traj_y_12);

                        // Store the CVT angles and corresponding edge values
                        theta_CVT_mc.push_back(mc_theta_CVT_value);
                        phi_CVT_mc.push_back(mc_phi_CVT_value);
                        theta_in_degrees_mc.push_back(**mc_theta);

                        // Store the corresponding edge values
                        traj_edge_1_mc.push_back(**mc_traj_edge_1);
                        traj_edge_3_mc.push_back(**mc_traj_edge_3);
                        traj_edge_5_mc.push_back(**mc_traj_edge_5);
                        traj_edge_7_mc.push_back(**mc_traj_edge_7);
                        traj_edge_12_mc.push_back(**mc_traj_edge_12);
                    }
                }
            }
        }

        // Create and fill histograms for theta_CVT vs theta and phi_CVT vs theta_CVT
        TH2D* h_theta_vs_theta_data_before = new TH2D("h_theta_vs_theta_data_before", ("#theta_{CVT} vs #theta Before Cuts (Data, " + particle_latex + ")").c_str(), nBins, 25, 150, nBins, 30, 150);
        h_theta_vs_theta_data_before->GetXaxis()->SetTitle("#theta");
        h_theta_vs_theta_data_before->GetYaxis()->SetTitle("#theta_{CVT}");

        TH2D* h_phi_vs_theta_CVT_data_before = new TH2D("h_phi_vs_theta_CVT_data_before", ("#phi_{CVT} vs #theta_{CVT} Before Cuts (Data, " + particle_latex + ")").c_str(), nBins, 0, 360, nBins, 25, 150);
        h_phi_vs_theta_CVT_data_before->GetXaxis()->SetTitle("#phi_{CVT}");
        h_phi_vs_theta_CVT_data_before->GetYaxis()->SetTitle("#theta_{CVT}");

        TH2D* h_theta_vs_theta_data_after = new TH2D("h_theta_vs_theta_data_after", ("#theta_{CVT} vs #theta After Cuts (Data, " + particle_latex + ")").c_str(), nBins, 25, 150, nBins, 30, 150);
        h_theta_vs_theta_data_after->GetXaxis()->SetTitle("#theta");
        h_theta_vs_theta_data_after->GetYaxis()->SetTitle("#theta_{CVT}");

        TH2D* h_phi_vs_theta_CVT_data_after = new TH2D("h_phi_vs_theta_CVT_data_after", ("#phi_{CVT} vs #theta_{CVT} After Cuts (Data, " + particle_latex + ")").c_str(), nBins, 0, 360, nBins, 25, 150);
        h_phi_vs_theta_CVT_data_after->GetXaxis()->SetTitle("#phi_{CVT}");
        h_phi_vs_theta_CVT_data_after->GetYaxis()->SetTitle("#theta_{CVT}");

        // Fill histograms using stored values for data
        for (size_t i = 0; i < theta_CVT_data.size(); ++i) {
            h_theta_vs_theta_data_before->Fill(theta_in_degrees_data[i], theta_CVT_data[i]);
            h_phi_vs_theta_CVT_data_before->Fill(phi_CVT_data[i], theta_CVT_data[i]);

            // Apply fiducial cut based on stored edge values
            if (cvt_fiducial(traj_edge_1_data[i], traj_edge_3_data[i], traj_edge_5_data[i], traj_edge_7_data[i], traj_edge_12_data[i])) {
                h_theta_vs_theta_data_after->Fill(theta_in_degrees_data[i], theta_CVT_data[i]);
                    h_phi_vs_theta_CVT_data_after->Fill(phi_CVT_data[i], theta_CVT_data[i]);
            }
        }

        TH2D* h_theta_vs_theta_mc_before = nullptr;
        TH2D* h_phi_vs_theta_CVT_mc_before = nullptr;
        TH2D* h_theta_vs_theta_mc_after = nullptr;
        TH2D* h_phi_vs_theta_CVT_mc_after = nullptr;

        // Create and fill histograms using stored values for MC
        if (mcReader) {
            h_theta_vs_theta_mc_before = new TH2D("h_theta_vs_theta_mc_before", ("#theta_{CVT} vs #theta Before Cuts (MC, " + particle_latex + ")").c_str(), nBins, 25, 150, nBins, 30, 150);
            h_theta_vs_theta_mc_before->GetXaxis()->SetTitle("#theta");
            h_theta_vs_theta_mc_before->GetYaxis()->SetTitle("#theta_{CVT}");

            h_phi_vs_theta_CVT_mc_before = new TH2D("h_phi_vs_theta_CVT_mc_before", ("#phi_{CVT} vs #theta_{CVT} Before Cuts (MC, " + particle_latex + ")").c_str(), nBins, 0, 360, nBins, 25, 150);
            h_phi_vs_theta_CVT_mc_before->GetXaxis()->SetTitle("#phi_{CVT}");
            h_phi_vs_theta_CVT_mc_before->GetYaxis()->SetTitle("#theta_{CVT}");

            h_theta_vs_theta_mc_after = new TH2D("h_theta_vs_theta_mc_after", ("#theta_{CVT} vs #theta After Cuts (MC, " + particle_latex + ")").c_str(), nBins, 25, 150, nBins, 30, 150);
            h_theta_vs_theta_mc_after->GetXaxis()->SetTitle("#theta");
            h_theta_vs_theta_mc_after->GetYaxis()->SetTitle("#theta_{CVT}");

            h_phi_vs_theta_CVT_mc_after = new TH2D("h_phi_vs_theta_CVT_mc_after", ("#phi_{CVT} vs #theta_{CVT} After Cuts (MC, " + particle_latex + ")").c_str(), nBins, 0, 360, nBins, 25, 150);
            h_phi_vs_theta_CVT_mc_after->GetXaxis()->SetTitle("#phi_{CVT}");
            h_phi_vs_theta_CVT_mc_after->GetYaxis()->SetTitle("#theta_{CVT}");

            for (size_t i = 0; i < theta_CVT_mc.size(); ++i) {
                h_theta_vs_theta_mc_before->Fill(theta_in_degrees_mc[i], theta_CVT_mc[i]);
                h_phi_vs_theta_CVT_mc_before->Fill(phi_CVT_mc[i], theta_CVT_mc[i]);

                // Apply fiducial cut based on stored edge values
                if (cvt_fiducial(traj_edge_1_mc[i], traj_edge_3_mc[i], traj_edge_5_mc[i], traj_edge_7_mc[i], traj_edge_12_mc[i])) {
                    h_theta_vs_theta_mc_after->Fill(theta_in_degrees_mc[i], theta_CVT_mc[i]);
                    h_phi_vs_theta_CVT_mc_after->Fill(phi_CVT_mc[i], theta_CVT_mc[i]);
                }
            }
        }

        // Create canvases for theta_CVT vs theta and phi_CVT vs theta_CVT
        TCanvas* c_theta_vs_theta_data = new TCanvas(("c_theta_vs_theta_data_" + particle_name).c_str(), ("#theta_{CVT} vs #theta and #phi_{CVT} vs #theta_{CVT} (Data, " + particle_latex + ")").c_str(), 1200, 1200);
        c_theta_vs_theta_data->Divide(2, 2);

        TCanvas* c_theta_vs_theta_mc = nullptr;
        if (mcReader) {
            c_theta_vs_theta_mc = new TCanvas(("c_theta_vs_theta_mc_" + particle_name).c_str(), ("#theta_{CVT} vs #theta and #phi_{CVT} vs #theta_{CVT} (MC, " + particle_latex + ")").c_str(), 1200, 1200);
            c_theta_vs_theta_mc->Divide(2, 2);
        }

        // Draw and save the theta_CVT vs theta and phi_CVT vs theta_CVT canvases (data and MC)
        c_theta_vs_theta_data->cd(1);
        gPad->SetMargin(0.2, 0.05, 0.1, 0.1);
        h_theta_vs_theta_data_before->Draw("COLZ");

        c_theta_vs_theta_data->cd(2);
        gPad->SetMargin(0.2, 0.05, 0.1, 0.1);
        h_phi_vs_theta_CVT_data_before->Draw("COLZ");

        c_theta_vs_theta_data->cd(3);
        gPad->SetMargin(0.2, 0.05, 0.1, 0.1);
        h_theta_vs_theta_data_after->Draw("COLZ");

        c_theta_vs_theta_data->cd(4);
        gPad->SetMargin(0.2, 0.05, 0.1, 0.1);
        h_phi_vs_theta_CVT_data_after->Draw("COLZ");

        c_theta_vs_theta_data->SaveAs(("output/calibration/cvt/positions/theta_vs_theta_data_" + particle_name + ".png").c_str());

        if (mcReader) {
            c_theta_vs_theta_mc->cd(1);
            gPad->SetMargin(0.2, 0.05, 0.1, 0.1);
            h_theta_vs_theta_mc_before->Draw("COLZ");

            c_theta_vs_theta_mc->cd(2);
            gPad->SetMargin(0.2, 0.05, 0.1, 0.1);
            h_phi_vs_theta_CVT_mc_before->Draw("COLZ");

            c_theta_vs_theta_mc->cd(3);
            gPad->SetMargin(0.2, 0.05, 0.1, 0.1);
            h_theta_vs_theta_mc_after->Draw("COLZ");

            c_theta_vs_theta_mc->cd(4);
            gPad->SetMargin(0.2, 0.05, 0.1, 0.1);
            h_phi_vs_theta_CVT_mc_after->Draw("COLZ");

            c_theta_vs_theta_mc->SaveAs(("output/calibration/cvt/positions/theta_vs_theta_mc_" + particle_name + ".png").c_str());
        }

        // Clean up
        delete h_theta_vs_theta_data_before;
        delete h_phi_vs_theta_CVT_data_before;
        delete h_theta_vs_theta_data_after;
        delete h_phi_vs_theta_CVT_data_after;

        if (mcReader) {
            delete h_theta_vs_theta_mc_before;
            delete h_phi_vs_theta_CVT_mc_before;
            delete h_theta_vs_theta_mc_after;
            delete h_phi_vs_theta_CVT_mc_after;
            delete c_theta_vs_theta_mc;
        }
        delete c_theta_vs_theta_data;
    }
    // Clean up the dynamically allocated memory for edge variables
    if (mc_traj_edge_1) delete mc_traj_edge_1;
    if (mc_traj_edge_3) delete mc_traj_edge_3;
    if (mc_traj_edge_5) delete mc_traj_edge_5;
    if (mc_traj_edge_7) delete mc_traj_edge_7;
    if (mc_traj_edge_12) delete mc_traj_edge_12;

    for (auto& ptr : mc_traj_x) {
        delete ptr;
    }
    for (auto& ptr : mc_traj_y) {
        delete ptr;
    }
}

// Helper function to check if a track is FD
bool is_fd_track(double track_sector_5) {
    return track_sector_5 != -9999;
}

// Helper function to check if a track is CD
bool is_cd_track(double track_sector_6) {
    return track_sector_6 != -9999;
}

// Helper function to fill and save histograms for each particle type
void fill_and_save_histograms(const std::map<int, std::pair<std::string, std::pair<TH2D*, TH2D*>>>& histograms, const std::string& dataset) {
    for (const auto& entry : histograms) {
        int pid = entry.first;
        const std::string& particle_name = entry.second.first;
        TH2D* h_fd = entry.second.second.first;
        TH2D* h_cd = entry.second.second.second;

        // Create a canvas with 1x2 subplots
        TCanvas* c = new TCanvas(("c_" + particle_name).c_str(), ("Energy Loss Distributions: " + dataset + ", " + particle_name).c_str(), 1600, 800);
        c->Divide(2, 1);

        // Add a title to the canvas
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.SetTextAlign(13);  // Align at top left
        latex.DrawLatexNDC(0.45, 0.97, (dataset + ", " + particle_name).c_str());

        // Plot FD
        c->cd(1);
        gPad->SetMargin(0.15, 0.15, 0.10, 0.1);  // Left, right, bottom, top margins
        gPad->SetLogz();
        h_fd->Draw("COLZ");

        // Plot CD
        c->cd(2);
        gPad->SetMargin(0.15, 0.15, 0.10, 0.1);  // Left, right, bottom, top margins
        gPad->SetLogz();
        h_cd->Draw("COLZ");

        // Save the canvas
        c->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/energy_loss_distributions_" + particle_name + ".png").c_str());

        // Clean up
        delete c;
        delete h_fd;
        delete h_cd;
    }
}

// Main energy loss distribution function
void energy_loss_distributions(TTreeReader& mcReader, const std::string& dataset) {
    // Particle types and their corresponding LaTeX names and x-axis ranges
    std::map<int, std::tuple<std::string, double, double>> particle_types = {
        // {11, {"e^{-}", 0.0, 7.0}},
        // {211, {"#pi^{+}", 0.0, 5.0}},
        // {-211, {"#pi^{-}", 0.0, 5.0}},
        // {321, {"k^{+}", 0.0, 5.0}},
        // {-321, {"k^{-}", 0.0, 5.0}},
        {2212, {"p", 0.0, 3.0}}
    };

    // Create histograms for each particle type
    std::map<int, std::pair<std::string, std::pair<TH2D*, TH2D*>>> histograms;
    for (const auto& particle : particle_types) {
        int pid = particle.first;
        const std::string& particle_name = std::get<0>(particle.second);
        double xMin = std::get<1>(particle.second);
        double xMax = std::get<2>(particle.second);

        TH2D* h_fd = new TH2D(("h_fd_" + particle_name).c_str(), ("(FD), " + particle_name).c_str(), 100, xMin, xMax, 100, -0.05, 0.10);
        TH2D* h_cd = new TH2D(("h_cd_" + particle_name).c_str(), ("(CD), " + particle_name).c_str(), 100, xMin, xMax, 100, -0.05, 0.10);

        h_fd->GetXaxis()->SetTitle("p (GeV)");
        h_fd->GetYaxis()->SetTitle("#Delta p (GeV)");
        h_fd->SetStats(false);

        h_cd->GetXaxis()->SetTitle("p (GeV)");
        h_cd->GetYaxis()->SetTitle("#Delta p (GeV)");
        h_cd->SetStats(false);

        histograms[pid] = std::make_pair(particle_name, std::make_pair(h_fd, h_cd));
    }

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Set up TTreeReaderValues for necessary branches
    TTreeReaderValue<double> mc_p(mcReader, "mc_p");
    TTreeReaderValue<double> p(mcReader, "p");
    TTreeReaderValue<int> pid(mcReader, "particle_pid");
    TTreeReaderValue<int> track_sector_5(mcReader, "track_sector_5");
    TTreeReaderValue<int> track_sector_6(mcReader, "track_sector_6");

    // Edge variables for FD and CD fiducial cuts
    TTreeReaderValue<double> edge_6(mcReader, "traj_edge_6");
    TTreeReaderValue<double> edge_18(mcReader, "traj_edge_18");
    TTreeReaderValue<double> edge_36(mcReader, "traj_edge_36");
    TTreeReaderValue<double> edge_1(mcReader, "traj_edge_1");
    TTreeReaderValue<double> edge_3(mcReader, "traj_edge_3");
    TTreeReaderValue<double> edge_5(mcReader, "traj_edge_5");
    TTreeReaderValue<double> edge_7(mcReader, "traj_edge_7");
    TTreeReaderValue<double> edge_12(mcReader, "traj_edge_12");

    // Loop over events
    while (mcReader.Next()) {
        double delta_p = *mc_p - *p;

        // Check if the current particle type is one of interest
        if (histograms.find(*pid) != histograms.end()) {
            // Check if the track is FD or CD and fill the appropriate histogram
            if (is_fd_track(*track_sector_6)) {
                if (dc_fiducial(*edge_6, *edge_18, *edge_36, *pid)) {
                    histograms[*pid].second.first->Fill(*p, delta_p);  // FD
                }
            } else if (is_cd_track(*track_sector_5)) {
                if (cvt_fiducial(*edge_1, *edge_3, *edge_5, *edge_7, *edge_12)) {
                    histograms[*pid].second.second->Fill(*p, delta_p);  // CD
                }
            }
        }
    }

    // Save the histograms
    fill_and_save_histograms(histograms, dataset);
}

bool is_above_deltap_curve(double p, double delta_p) {
    return (delta_p > 0.011 / pow(p, 1.05));
}

bool is_above_theta_dc_curve(double p, double theta_dc_1) {
    double curve_value = -53.1468 + 79.6131 * pow(p - 0.3, 0.05739);
    return (theta_dc_1 > curve_value);
}

// Main FD-specific function
void energy_loss_fd_distributions(TTreeReader& mcReader, const std::string& dataset) {
    // Particle types and their corresponding LaTeX names and x-axis ranges
    std::map<int, std::tuple<std::string, double, double>> particle_types = {
        // {11, {"e^{-}", 0.0, 7.0}},
        // {211, {"#pi^{+}", 0.0, 5.0}},
        // {-211, {"#pi^{-}", 0.0, 5.0}},
        // {321, {"k^{+}", 0.0, 5.0}},
        // {-321, {"k^{-}", 0.0, 5.0}},
        {2212, {"p", 0.0, 4.0}}
    };

    // Create histograms for each particle type
    std::map<int, std::vector<TH2D*>> histograms;
    for (const auto& particle : particle_types) {
        int pid = particle.first;
        const std::string& particle_name = std::get<0>(particle.second);
        double xMin = std::get<1>(particle.second);
        double xMax = std::get<2>(particle.second);

        histograms[pid] = {
            new TH2D(("h_above_deltap_" + particle_name).c_str(), ("(Above), " + particle_name).c_str(), 75, xMin, xMax, 75, -0.05, 0.10),
            new TH2D(("h_above_thetadc1_" + particle_name).c_str(), ("(Above), " + particle_name).c_str(), 75, xMin, xMax, 75, 0, 40),  // Adjusted to 0-40 degrees
            new TH2D(("h_above_theta_" + particle_name).c_str(), ("(Above), " + particle_name).c_str(), 75, xMin, xMax, 75, 0, 40),  // Adjusted to 0-40 degrees
            new TH2D(("h_below_deltap_" + particle_name).c_str(), ("(Below), " + particle_name).c_str(), 75, xMin, xMax, 75, -0.05, 0.10),
            new TH2D(("h_below_thetadc1_" + particle_name).c_str(), ("(Below), " + particle_name).c_str(), 75, xMin, xMax, 75, 0, 40),  // Adjusted to 0-40 degrees
            new TH2D(("h_below_theta_" + particle_name).c_str(), ("(Below), " + particle_name).c_str(), 75, xMin, xMax, 75, 0, 40)  // Adjusted to 0-40 degrees
        };

        // Set axis labels
        histograms[pid][0]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][0]->GetYaxis()->SetTitle("#Delta p (GeV)");
        histograms[pid][1]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][1]->GetYaxis()->SetTitle("#theta_{DC_{region 1}} (degrees)");
        histograms[pid][2]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][2]->GetYaxis()->SetTitle("#theta (degrees)");
        histograms[pid][3]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][3]->GetYaxis()->SetTitle("#Delta p (GeV)");
        histograms[pid][4]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][4]->GetYaxis()->SetTitle("#theta_{DC_{region 1}} (degrees)");
        histograms[pid][5]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][5]->GetYaxis()->SetTitle("#theta (degrees)");

        for (auto& hist : histograms[pid]) {
            hist->SetStats(false);
        }
    }

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Set up TTreeReaderValues for necessary branches
    TTreeReaderValue<double> mc_p(mcReader, "mc_p");
    TTreeReaderValue<double> p(mcReader, "p");
    TTreeReaderValue<double> traj_x_6(mcReader, "traj_x_6");
    TTreeReaderValue<double> traj_y_6(mcReader, "traj_y_6");
    TTreeReaderValue<double> traj_z_6(mcReader, "traj_z_6");
    TTreeReaderValue<double> theta(mcReader, "theta");
    TTreeReaderValue<int> pid(mcReader, "particle_pid");
    TTreeReaderValue<int> track_sector_6(mcReader, "track_sector_6");

    // Edge variables for FD fiducial cuts
    TTreeReaderValue<double> edge_6(mcReader, "traj_edge_6");
    TTreeReaderValue<double> edge_18(mcReader, "traj_edge_18");
    TTreeReaderValue<double> edge_36(mcReader, "traj_edge_36");

    // Loop over events
    // for (int i = 0; i < 1e7; ++i) {
    //     mcReader.Next();
    while (mcReader.Next()) {
        double delta_p = *mc_p - *p;
        double theta_dc_1 = calculate_theta(*traj_x_6, *traj_y_6, *traj_z_6);

        // Check if the current particle type is one of interest
        if (histograms.find(*pid) != histograms.end()) { 
            bool above_curve = is_above_deltap_curve(*p, delta_p);

            if (is_fd_track(*track_sector_6)) {
                if (dc_fiducial(*edge_6, *edge_18, *edge_36, *pid)) {
                    int index_deltap = above_curve ? 0 : 3;
                    int index_thetadc1 = above_curve ? 1 : 4;
                    int index_theta = above_curve ? 2 : 5;

                    // Fill histograms
                    histograms[*pid][index_deltap]->Fill(*p, delta_p);
                    histograms[*pid][index_thetadc1]->Fill(*p, theta_dc_1);
                    histograms[*pid][index_theta]->Fill(*p, *theta);
                }
            }
        }
    }

    for (const auto& entry : histograms) {
        int pid = entry.first;
        const std::string& particle_name = std::get<0>(particle_types[pid]);

        // Create a canvas with 2x3 subplots
        TCanvas* c = new TCanvas(("c_fd_" + particle_name).c_str(), ("FD Energy Loss Distributions: " + dataset + ", " + particle_name).c_str(), 1800, 1200);
        c->Divide(3, 2);

        // Move the LaTeX title up to avoid clipping
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.SetTextAlign(11);  // Align at top left
        latex.DrawLatexNDC(0.42, 0.525, (dataset + ", " + particle_name).c_str());

        // Define the curve function
        TF1* pass1_curve = new TF1("pass1 curve", "0.088/pow(x, 1.5)", 0.1, std::get<2>(particle_types[pid]));
        pass1_curve->SetLineColor(kRed);
        pass1_curve->SetLineWidth(4);
        TF1* pass2_curve = new TF1("pass2 curve", "0.011/pow(x, 1.05)", 0.1, std::get<2>(particle_types[pid]));
        pass2_curve->SetLineColor(kBlack);
        pass2_curve->SetLineWidth(4);

        // Define the new curve based on the provided formula
        TF1* new_curve_region1 = new TF1("new_curve_region1", "-53.1468 + 79.6131*pow(x-0.3,0.05739)", 0.30, std::get<2>(particle_types[pid]));
        new_curve_region1->SetLineColor(kRed);
        new_curve_region1->SetLineWidth(4);
        new_curve_region1->SetNpx(5000);  // Increase the number of points along the curve

        // Plot histograms and curve in the correct order
        c->cd(1);
        gPad->SetMargin(0.15, 0.15, 0.20, 0.1);  // Left, right, bottom, top margins
        gPad->SetLogz();
        entry.second[0]->Draw("COLZ");
        pass1_curve->Draw("same");
        pass2_curve->Draw("same");

        c->cd(2);
        gPad->SetMargin(0.15, 0.15, 0.20, 0.1);
        gPad->SetLogz();
        entry.second[1]->Draw("COLZ");
        new_curve_region1->Draw("same");

        c->cd(3);
        gPad->SetMargin(0.15, 0.15, 0.20, 0.1);
        gPad->SetLogz();
        entry.second[2]->Draw("COLZ");

        c->cd(4);
        gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
        gPad->SetLogz();
        entry.second[3]->Draw("COLZ");
        pass1_curve->Draw("same");
        pass2_curve->Draw("same");

        c->cd(5);
        gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
        gPad->SetLogz();
        entry.second[4]->Draw("COLZ");
        new_curve_region1->Draw("same");

        c->cd(6);
        gPad->SetMargin(0.15, 0.15 , 0.1, 0.1);
        gPad->SetLogz();
        entry.second[5]->Draw("COLZ");

        // Save the canvas
        c->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/fd_energy_loss_distributions_" + particle_name + ".png").c_str());

        // Clean up
        delete c;
        delete pass1_curve;
        delete pass2_curve;
        for (auto& hist : entry.second) {
            delete hist;
        }
    }
}

void energy_loss_fd_distributions_theta_dc(TTreeReader& mcReader, const std::string& dataset) {
    // Particle types and their corresponding LaTeX names and x-axis ranges
    std::map<int, std::tuple<std::string, double, double>> particle_types = {
        // {11, {"e^{-}", 0.0, 7.0}},
        // {211, {"#pi^{+}", 0.0, 5.0}},
        // {-211, {"#pi^{-}", 0.0, 5.0}},
        // {321, {"k^{+}", 0.0, 5.0}},
        // {-321, {"k^{-}", 0.0, 5.0}},
        {2212, {"p", 0.0, 4.0}}
    };

    // Create histograms for each particle type
    std::map<int, std::vector<TH2D*>> histograms;
    for (const auto& particle : particle_types) {
        int pid = particle.first;
        const std::string& particle_name = std::get<0>(particle.second);
        double xMin = std::get<1>(particle.second);
        double xMax = std::get<2>(particle.second);

        histograms[pid] = {
            new TH2D(("h_above_deltap_" + particle_name).c_str(), ("(Above), " + particle_name).c_str(), 75, xMin, xMax, 75, -0.05, 0.10),
            new TH2D(("h_above_thetadc1_" + particle_name).c_str(), ("(Above), " + particle_name).c_str(), 75, xMin, xMax, 75, 0, 40),  
            new TH2D(("h_above_theta_" + particle_name).c_str(), ("(Above), " + particle_name).c_str(), 75, xMin, xMax, 75, 0, 40),  
            new TH2D(("h_below_deltap_" + particle_name).c_str(), ("(Below), " + particle_name).c_str(), 75, xMin, xMax, 75, -0.05, 0.10),
            new TH2D(("h_below_thetadc1_" + particle_name).c_str(), ("(Below), " + particle_name).c_str(), 75, xMin, xMax, 75, 0, 40),  
            new TH2D(("h_below_theta_" + particle_name).c_str(), ("(Below), " + particle_name).c_str(), 75, xMin, xMax, 75, 0, 40)  
        };

        // Set axis labels
        histograms[pid][0]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][0]->GetYaxis()->SetTitle("#Delta p (GeV)");
        histograms[pid][1]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][1]->GetYaxis()->SetTitle("#theta_{DC_{region 1}} (degrees)");
        histograms[pid][2]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][2]->GetYaxis()->SetTitle("#theta (degrees)");
        histograms[pid][3]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][3]->GetYaxis()->SetTitle("#Delta p (GeV)");
        histograms[pid][4]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][4]->GetYaxis()->SetTitle("#theta_{DC_{region 1}} (degrees)");
        histograms[pid][5]->GetXaxis()->SetTitle("p (GeV)"); histograms[pid][5]->GetYaxis()->SetTitle("#theta (degrees)");

        for (auto& hist : histograms[pid]) {
            hist->SetStats(false);
        }
    }

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Set up TTreeReaderValues for necessary branches
    TTreeReaderValue<double> mc_p(mcReader, "mc_p");
    TTreeReaderValue<double> p(mcReader, "p");
    TTreeReaderValue<double> traj_x_6(mcReader, "traj_x_6");
    TTreeReaderValue<double> traj_y_6(mcReader, "traj_y_6");
    TTreeReaderValue<double> traj_z_6(mcReader, "traj_z_6");
    TTreeReaderValue<double> theta(mcReader, "theta");
    TTreeReaderValue<int> pid(mcReader, "particle_pid");
    TTreeReaderValue<int> track_sector_6(mcReader, "track_sector_6");

    // Edge variables for FD fiducial cuts
    TTreeReaderValue<double> edge_6(mcReader, "traj_edge_6");
    TTreeReaderValue<double> edge_18(mcReader, "traj_edge_18");
    TTreeReaderValue<double> edge_36(mcReader, "traj_edge_36");

    // Loop over events
    while (mcReader.Next()) {
        double delta_p = *mc_p - *p;
        double theta_dc_1 = calculate_theta(*traj_x_6, *traj_y_6, *traj_z_6);

        // Check if the current particle type is one of interest
        if (histograms.find(*pid) != histograms.end()) { 
            bool above_curve = is_above_theta_dc_curve(*p, theta_dc_1);

            if (is_fd_track(*track_sector_6)) {
                if (dc_fiducial(*edge_6, *edge_18, *edge_36, *pid)) {
                    int index_deltap = above_curve ? 0 : 3;
                    int index_thetadc1 = above_curve ? 1 : 4;
                    int index_theta = above_curve ? 2 : 5;

                    // Fill histograms
                    histograms[*pid][index_deltap]->Fill(*p, delta_p);
                    histograms[*pid][index_thetadc1]->Fill(*p, theta_dc_1);
                    histograms[*pid][index_theta]->Fill(*p, *theta);
                }
            }
        }
    }

    for (const auto& entry : histograms) {
        int pid = entry.first;
        const std::string& particle_name = std::get<0>(particle_types[pid]);

        // Create a canvas with 2x3 subplots
        TCanvas* c = new TCanvas(("c_fd_" + particle_name).c_str(), ("FD Energy Loss Distributions (theta DC cut): " + dataset + ", " + particle_name).c_str(), 1800, 1200);
        c->Divide(3, 2);

        // Move the LaTeX title up to avoid clipping
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.SetTextAlign(11);  // Align at top left
        latex.DrawLatexNDC(0.42, 0.525, (dataset + ", " + particle_name).c_str());

        // Define the curve function
        TF1* pass1_curve = new TF1("pass1 curve", "0.088/pow(x, 1.5)", 0.1, std::get<2>(particle_types[pid]));
        pass1_curve->SetLineColor(kRed);
        pass1_curve->SetLineWidth(4);
        TF1* pass2_curve = new TF1("pass2 curve", "0.011/pow(x, 1.05)", 0.1, std::get<2>(particle_types[pid]));
        pass2_curve->SetLineColor(kBlack);
        pass2_curve->SetLineWidth(4);

        // Define the new curve based on the provided formula
        TF1* new_curve_region1 = new TF1("new_curve_region1", "-53.1468 + 79.6131*pow(x-0.3,0.05739)", 0.30, std::get<2>(particle_types[pid]));
        new_curve_region1->SetLineColor(kRed);
        new_curve_region1->SetLineWidth(4);
        new_curve_region1->SetNpx(5000);  // Increase the number of points along the curve

        // Plot histograms and curve in the correct order
        c->cd(1);
        gPad->SetMargin(0.15, 0.15, 0.20, 0.1);  // Left, right, bottom, top margins
        gPad->SetLogz();
        entry.second[0]->Draw("COLZ");
        pass1_curve->Draw("same");
        pass2_curve->Draw("same");

        c->cd(2);
        gPad->SetMargin(0.15, 0.15, 0.20, 0.1);
        gPad->SetLogz();
        entry.second[1]->Draw("COLZ");
        new_curve_region1->Draw("same");
        c->cd(3);
        gPad->SetMargin(0.15, 0.15, 0.20, 0.1);
        gPad->SetLogz();
        entry.second[2]->Draw("COLZ");

        c->cd(4);
        gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
        gPad->SetLogz();
        entry.second[3]->Draw("COLZ");
        pass1_curve->Draw("same");
        pass2_curve->Draw("same");

        c->cd(5);
        gPad->SetMargin(0.15, 0.15, 0.1, 0.1);
        gPad->SetLogz();
        entry.second[4]->Draw("COLZ");
        new_curve_region1->Draw("same");

        c->cd(6);
        gPad->SetMargin(0.15, 0.15 , 0.1, 0.1);
        gPad->SetLogz();
        entry.second[5]->Draw("COLZ");

        // Save the canvas
        c->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/fd_energy_loss_distributions_theta_dc_" + particle_name + ".png").c_str());

        // Clean up
        delete c;
        delete pass1_curve;
        delete pass2_curve;
        for (auto& hist : entry.second) {
            delete hist;
        }
    }
}

void plot_and_fit_parameters(const std::vector<std::pair<double, double>>& theta_bins,
                             const std::vector<double>& A_values,
                             const std::vector<double>& A_errors,
                             const std::vector<double>& B_values,
                             const std::vector<double>& B_errors,
                             const std::vector<double>& C_values,
                             const std::vector<double>& C_errors,
                             const std::string& particle_name,
                             const std::string& dataset,
                             const std::string& prefix) {
    // Create a new canvas for the fitted parameters
    TCanvas* c_fit_params = new TCanvas(("c_fit_params_" + particle_name).c_str(), 
                                        ("Fit Parameters: " + dataset + ", " + particle_name).c_str(), 
                                        1600, 800);
    c_fit_params->Divide(3, 1);  // 1 row, 2 columns

    // Plot A(#theta)
    c_fit_params->cd(1);
    TGraphErrors* graph_A = new TGraphErrors(theta_bins.size());
    for (size_t i = 0; i < theta_bins.size(); ++i) {
        double theta_midpoint = 0.5 * (theta_bins[i].first + theta_bins[i].second);
        graph_A->SetPoint(i, theta_midpoint, A_values[i]);
        graph_A->SetPointError(i, 0.0, A_errors[i]);
    }
    
    if (prefix == "p") {
        graph_A->GetYaxis()->SetRangeUser(-0.02, 0.02);  // Set y-axis range
        graph_A->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_A->SetTitle(("A_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", FD;#theta (degrees);A_{" + prefix + "}(#theta) (GeV)").c_str());
    } else {
        graph_A->GetYaxis()->SetRangeUser(-0.5, 0.5);  // Set y-axis range
        graph_A->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_A->SetTitle(("A_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", FD;#theta (degrees);A_{" + prefix + "}(#theta)").c_str());
    }
    
    graph_A->SetMarkerStyle(20);  // Set marker style to a filled circle
    gPad->SetLeftMargin(0.2);  // Increase left margin
    graph_A->Draw("AP");

    // Fit A(#theta) to a 2nd order polynomial
    TF1* fit_A = new TF1("fit_A", "[0]+[1]*x+[2]*x*x", theta_bins.front().first, theta_bins.back().second);
    graph_A->Fit(fit_A, "Q");  // Silent fit
    fit_A->Draw("same");

    // Add fit results and chi2/ndf to the plot
    TPaveText* pt_A = new TPaveText(0.7, 0.75, 0.9, 0.9, "NDC");
    pt_A->AddText(Form("p0 = %.7f", fit_A->GetParameter(0)));
    pt_A->AddText(Form("p1 = %.7f", fit_A->GetParameter(1)));
    pt_A->AddText(Form("p2 = %.7f", fit_A->GetParameter(2)));
    pt_A->AddText(Form("#chi^{2}/ndf = %.3f", fit_A->GetChisquare() / fit_A->GetNDF()));
    pt_A->SetBorderSize(1);
    pt_A->SetFillColor(0);
    pt_A->Draw();

    // Plot B(#theta)
    c_fit_params->cd(2);
    TGraphErrors* graph_B = new TGraphErrors(theta_bins.size());
    for (size_t i = 0; i < theta_bins.size(); ++i) {
        double theta_midpoint = 0.5 * (theta_bins[i].first + theta_bins[i].second);
        graph_B->SetPoint(i, theta_midpoint, B_values[i]);
        graph_B->SetPointError(i, 0.0, B_errors[i]);
    }
    
    if (prefix == "p") {
        graph_B->GetYaxis()->SetRangeUser(-0.02, 0.02);  // Set y-axis range
        graph_B->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_B->SetTitle(("B_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", FD;#theta (degrees);B_{" + prefix + "}(#theta) (GeV^{2})").c_str());
    } else {
        graph_B->GetYaxis()->SetRangeUser(-0.5, 0.5);  // Set y-axis range
        graph_B->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_B->SetTitle(("B_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", FD;#theta (degrees);B_{" + prefix + "}(#theta)").c_str());
    }
    graph_B->SetMarkerStyle(20);  // Set marker style to a filled circle
    gPad->SetLeftMargin(0.2);  // Increase left margin
    graph_B->Draw("AP");

    // Fit B(#theta) to a 2nd order polynomial
    TF1* fit_B = new TF1("fit_B", "[0]+[1]*x+[2]*x*x", theta_bins.front().first, theta_bins.back().second);
    graph_B->Fit(fit_B, "Q");  // Silent fit
    fit_B->Draw("same");

    // Add fit results and chi2/ndf to the plot
    TPaveText* pt_B = new TPaveText(0.7, 0.75, 0.9, 0.9, "NDC");
    pt_B->AddText(Form("p0 = %.7f", fit_B->GetParameter(0)));
    pt_B->AddText(Form("p1 = %.7f", fit_B->GetParameter(1)));
    pt_B->AddText(Form("p2 = %.7f", fit_B->GetParameter(2)));
    pt_B->AddText(Form("#chi^{2}/ndf = %.3f", fit_B->GetChisquare() / fit_B->GetNDF()));
    pt_B->SetBorderSize(1);
    pt_B->SetFillColor(0);
    pt_B->Draw();

    // Plot C(#theta)
    c_fit_params->cd(3);
    TGraphErrors* graph_C = new TGraphErrors(theta_bins.size());
    for (size_t i = 0; i < theta_bins.size(); ++i) {
        double theta_midpoint = 0.5 * (theta_bins[i].first + theta_bins[i].second);
        graph_C->SetPoint(i, theta_midpoint, C_values[i]);
        graph_C->SetPointError(i, 0.0, C_errors[i]);
    }
    if (prefix == "p") {
        graph_C->GetYaxis()->SetRangeUser(-0.02, 0.02);  // Set y-axis range
        graph_C->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_C->SetTitle(("C_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", FD;#theta (degrees);C_{" + prefix + "}(#theta) (GeV^{3})").c_str());
    } else {
        graph_C->GetYaxis()->SetRangeUser(-0.5, 0.5);  // Set y-axis range
        graph_C->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_C->SetTitle(("C_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", FD;#theta (degrees);C_{" + prefix + "}(#theta)").c_str());
    }
    graph_C->SetMarkerStyle(20);  // Set marker style to a filled circle
    gPad->SetLeftMargin(0.2);  // Increase left margin
    graph_C->Draw("AP");

    // Fit C(#theta) to a 2nd order polynomial
    TF1* fit_C = new TF1("fit_C", "[0]+[1]*x+[2]*x*x", theta_bins.front().first, theta_bins.back().second);
    graph_C->Fit(fit_C, "Q");  // Silent fit
    fit_C->Draw("same");

    // Add fit results and chi2/ndf to the plot
    TPaveText* pt_C = new TPaveText(0.7, 0.75, 0.9, 0.9, "NDC");
    pt_C->AddText(Form("p0 = %.7f", fit_C->GetParameter(0)));
    pt_C->AddText(Form("p1 = %.7f", fit_C->GetParameter(1)));
    pt_C->AddText(Form("p2 = %.7f", fit_C->GetParameter(2)));
    pt_C->AddText(Form("#chi^{2}/ndf = %.3f", fit_C->GetChisquare() / fit_C->GetNDF()));
    pt_C->SetBorderSize(1);
    pt_C->SetFillColor(0);
    pt_C->Draw();

    // Print out the functional form of A(theta) in LaTeX format
    std::cout << "A_" << prefix << "(\\theta) = ";
    for (int i = 0; i <= 2; ++i) {
        double coeff = fit_A->GetParameter(i);
        if (i == 0) {
            std::cout << Form("%.7f", coeff);
        } else if (i == 1) {
            std::cout << Form(" %+.7f\\theta", coeff);
        } else {
            std::cout << Form(" %+.7f\\theta^%d", coeff, i);
        }
    }
    std::cout << std::endl;

    // Print out the functional form of B(theta) in LaTeX format
    std::cout << "B_" << prefix << "(\\theta) = ";
    for (int i = 0; i <= 2; ++i) {
        double coeff = fit_B->GetParameter(i);
        if (i == 0) {
            std::cout << Form("%.8f", coeff);
        } else if (i == 1) {
            std::cout << Form(" %+.8f\\theta", coeff);
        } else {
            std::cout << Form(" %+.8f\\theta^%d", coeff, i);
        }
    }
    std::cout << std::endl;

    // Print out the functional form of C(theta) in LaTeX format
    std::cout << "C_" << prefix << "(\\theta) = ";
    for (int i = 0; i <= 2; ++i) {
        double coeff = fit_C->GetParameter(i);
        if (i == 0) {
            std::cout << Form("%.8f", coeff);
        } else if (i == 1) {
            std::cout << Form(" %+.8f\\theta", coeff);
        } else {
            std::cout << Form(" %+.8f\\theta^%d", coeff, i);
        }
    }
    std::cout << std::endl;

    // Save the fit parameters canvas
    c_fit_params->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/fit_params_" + prefix + "_" + particle_name + ".png").c_str());

    // Clean up memory
    delete graph_A;
    delete graph_B;
    delete graph_C;
    delete fit_A;
    delete fit_B;
    delete fit_C;
    delete pt_A;
    delete pt_B;
    delete pt_C;
    delete c_fit_params;
}

void energy_loss_distributions_delta_p_fd(TTreeReader& mcReader, const std::string& dataset) {
    // Particle types and their corresponding LaTeX names and x-axis ranges
    std::map<int, std::tuple<std::string, double, double>> particle_types = {
        {2212, {"p", 0.0, 5.0}}
    };

    // // Define theta bins
    std::vector<std::pair<double, double>> theta_bins;
    if (dataset == "rga_fa18_out") {
        // Outbending case: 24 equally spaced bins between 8 and 44 degrees
        theta_bins = {
            {8.0, 9.5}, {9.5, 11.0}, {11.0, 12.5}, {12.5, 14.0},
            {14.0, 15.5}, {15.5, 17.0}, {17.0, 18.5}, {18.5, 20.0},
            {20.0, 21.5}, {21.5, 23.0}, {23.0, 24.5}, {24.5, 26.0},
            {26.0, 27.5}, {27.5, 29.0}, {29.0, 30.5}, {30.5, 32.0},
            {32.0, 33.5}, {33.5, 35.0}, {35.0, 36.5}, {36.5, 38.0},
            {38.0, 39.5}, {39.5, 41.0}, {41.0, 42.5}, {42.5, 44.0}
        };
    } else {
        // Inbending case (default): existing bins
        theta_bins = {
            {5.0, 6.3043}, {6.3043, 7.6087}, {7.6087, 8.9130},
            {8.9130, 10.2174}, {10.2174, 11.5217}, {11.5217, 12.8261},
            {12.8261, 14.1304}, {14.1304, 15.4348}, {15.4348, 16.7391},
            {16.7391, 18.0435}, {18.0435, 19.3478}, {19.3478, 20.6522},
            {20.6522, 21.9565}, {21.9565, 23.2609}, {23.2609, 24.5652},
            {24.5652, 25.8696}, {25.8696, 27.1739}, {27.1739, 28.4783},
            {28.4783, 29.7826}, {29.7826, 31.0870}, {31.0870, 32.3913},
            {32.3913, 33.6957}, {33.6957, 35.0}, {35.0, 42.0}
        };
    }

    // Create histograms for each particle type and theta bin
    std::map<int, std::vector<TH2D*>> histograms;
    for (const auto& particle : particle_types) {
        int pid = particle.first;
        const std::string& particle_name = std::get<0>(particle.second);
        double xMin = std::get<1>(particle.second);
        double xMax = std::get<2>(particle.second);

        histograms[pid].resize(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            std::string bin_label = TString::Format("#theta [%.1f, %.1f]", theta_bins[i].first, theta_bins[i].second).Data();

            histograms[pid][i] = new TH2D(
                ("h_deltap_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label.c_str(),
                75, xMin, xMax, 75, -0.05, 0.05);

            // Set axis labels
            histograms[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms[pid][i]->GetYaxis()->SetTitle("#Deltap");

            histograms[pid][i]->SetStats(false);
            histograms[pid][i]->GetXaxis()->SetLabelSize(0.04); // Increase font size for axes labels
            histograms[pid][i]->GetYaxis()->SetLabelSize(0.04);
        }
    }

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Set up TTreeReaderValues for necessary branches
    TTreeReaderValue<int> track_sector_6(mcReader, "track_sector_6");
    TTreeReaderValue<double> mc_p(mcReader, "mc_p");
    TTreeReaderValue<double> p(mcReader, "p");
    TTreeReaderValue<double> mc_theta(mcReader, "mc_theta");
    TTreeReaderValue<double> theta(mcReader, "theta");
    TTreeReaderValue<double> traj_x_6(mcReader, "traj_x_6");
    TTreeReaderValue<double> traj_y_6(mcReader, "traj_y_6");
    TTreeReaderValue<double> traj_z_6(mcReader, "traj_z_6");
    TTreeReaderValue<int> pid(mcReader, "particle_pid");

    // Edge variables for FD fiducial cuts
    TTreeReaderValue<double> edge_6(mcReader, "traj_edge_6");
    TTreeReaderValue<double> edge_18(mcReader, "traj_edge_18");
    TTreeReaderValue<double> edge_36(mcReader, "traj_edge_36");

    // Loop over events
    // for (int i = 0; i < 1e7; ++i) {
    //     mcReader.Next();
    while (mcReader.Next()) {
        if (!is_fd_track(*track_sector_6)) continue;
        if (!dc_fiducial(*edge_6, *edge_18, *edge_36, *pid)) continue;
        double delta_p = *mc_p - *p;
        double theta_dc_1 = calculate_theta(*traj_x_6, *traj_y_6, *traj_z_6);

        // Check if the current particle type is one of interest and if the track is below the curve
        // if (histograms.find(*pid) != histograms.end() && !is_above_theta_dc_curve(*p, theta_dc_1)) {
        if (histograms.find(*pid) != histograms.end() ) {
            for (size_t i = 0; i < theta_bins.size(); ++i) {
                if (*theta >= theta_bins[i].first && *theta < theta_bins[i].second) {
                    histograms[*pid][i]->Fill(*p, delta_p);
                    break;
                }
            }
        }
    }

    // Save the histograms into a canvas
    for (const auto& entry : histograms) {
        int pid = entry.first;
        const std::string& particle_name = std::get<0>(particle_types[pid]);

        TCanvas* c_deltap = new TCanvas(("c_deltap_" + particle_name).c_str(), ("Delta p Distributions: " + dataset + ", " + particle_name).c_str(), 2000, 1200);
        c_deltap->Divide(6, 4);  // 20 subplots

        std::vector<TF1*> fit_deltap(theta_bins.size());
        std::vector<double> A_values(theta_bins.size());
        std::vector<double> A_errors(theta_bins.size());
        std::vector<double> B_values(theta_bins.size());
        std::vector<double> B_errors(theta_bins.size());
        std::vector<double> C_values(theta_bins.size());
        std::vector<double> C_errors(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            // Ensure we are drawing on the correct pad
            c_deltap->cd(i + 1);
            gPad->SetMargin(0.15, 0.15, 0.20, 0.1);  // Left, right, bottom, top margins
            gPad->SetLogz();

            // Create profile histograms
            TProfile* prof_deltap = histograms[pid][i]->ProfileX();

            // Set the range of the profile to start from 0.5
            prof_deltap->GetXaxis()->SetRangeUser(0.5, std::get<2>(particle_types[pid]));

            // Fit the profiles with appropriate functions
            fit_deltap[i] = new TF1(("fit_deltap_" + std::to_string(i)).c_str(), "[0] + [1]/x + [2]/x^2", 0.5, std::get<2>(particle_types[pid]));
            
            prof_deltap->Fit(fit_deltap[i], "Q"); // Silent fit

            // Set the range of the fit function for plotting
            fit_deltap[i]->SetRange(0.5, std::get<2>(particle_types[pid]));

            // Store the fit parameters
            A_values[i] = fit_deltap[i]->GetParameter(0);
            A_errors[i] = fit_deltap[i]->GetParError(0);
            B_values[i] = fit_deltap[i]->GetParameter(1);
            B_errors[i] = fit_deltap[i]->GetParError(1);
            C_values[i] = fit_deltap[i]->GetParameter(2);
            C_errors[i] = fit_deltap[i]->GetParError(2);

            histograms[pid][i]->Draw("COLZ");
            prof_deltap->Draw("same");  // Draw the fit on top of the profile
        }

        // Add centered text "dataset, FD" on the canvas
        c_deltap->cd();  // Switch to the main canvas (not any specific pad)
        TLatex latex;
        latex.SetNDC();  // Use normalized coordinates (0,0) to (1,1)
        latex.SetTextSize(0.035);  // Set the text size
        latex.DrawLatex(0.425, 0.5, (dataset + ", FD").c_str());  // Add text in the center

        // Save the canvas
        c_deltap->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/delta_p_distributions_" + particle_name + ".png").c_str());

        // Use the new modular function for the fitted parameters
        plot_and_fit_parameters(theta_bins, A_values, A_errors, B_values, B_errors, C_values, C_errors, particle_name, dataset, "p");

        // Clean up memory
        for (size_t i = 0; i < theta_bins.size(); ++i) {
            delete fit_deltap[i];
            delete histograms[pid][i];
        }
        delete c_deltap;
    }
}

void energy_loss_distributions_delta_theta_fd(TTreeReader& mcReader, const std::string& dataset) {
    // Particle types and their corresponding LaTeX names and x-axis ranges
    std::map<int, std::tuple<std::string, double, double>> particle_types = {
        {2212, {"p", 0.0, 5.0}}
    };

    // // Define theta bins
    std::vector<std::pair<double, double>> theta_bins;
    if (dataset == "rga_fa18_out") {
        // Outbending case: 24 equally spaced bins between 8 and 44 degrees
        theta_bins = {
            {8.0, 9.5}, {9.5, 11.0}, {11.0, 12.5}, {12.5, 14.0},
            {14.0, 15.5}, {15.5, 17.0}, {17.0, 18.5}, {18.5, 20.0},
            {20.0, 21.5}, {21.5, 23.0}, {23.0, 24.5}, {24.5, 26.0},
            {26.0, 27.5}, {27.5, 29.0}, {29.0, 30.5}, {30.5, 32.0},
            {32.0, 33.5}, {33.5, 35.0}, {35.0, 36.5}, {36.5, 38.0},
            {38.0, 39.5}, {39.5, 41.0}, {41.0, 42.5}, {42.5, 44.0}
        };
    } else {
        // Inbending case (default): existing bins
        theta_bins = {
            {5.0, 6.3043}, {6.3043, 7.6087}, {7.6087, 8.9130},
            {8.9130, 10.2174}, {10.2174, 11.5217}, {11.5217, 12.8261},
            {12.8261, 14.1304}, {14.1304, 15.4348}, {15.4348, 16.7391},
            {16.7391, 18.0435}, {18.0435, 19.3478}, {19.3478, 20.6522},
            {20.6522, 21.9565}, {21.9565, 23.2609}, {23.2609, 24.5652},
            {24.5652, 25.8696}, {25.8696, 27.1739}, {27.1739, 28.4783},
            {28.4783, 29.7826}, {29.7826, 31.0870}, {31.0870, 32.3913},
            {32.3913, 33.6957}, {33.6957, 35.0}, {35.0, 42.0}
        };
    }

    // Create histograms for each particle type and theta bin
    std::map<int, std::vector<TH2D*>> histograms;
    for (const auto& particle : particle_types) {
        int pid = particle.first;
        const std::string& particle_name = std::get<0>(particle.second);
        double xMin = std::get<1>(particle.second);
        double xMax = std::get<2>(particle.second);

        histograms[pid].resize(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            std::string bin_label = TString::Format("#theta [%.1f, %.1f]", theta_bins[i].first, theta_bins[i].second).Data();

            histograms[pid][i] = new TH2D(
                ("h_deltatheta_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label.c_str(),
                75, xMin, xMax, 75, -1, 1);

            // Set axis labels
            histograms[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms[pid][i]->GetYaxis()->SetTitle("#Delta#theta");

            histograms[pid][i]->SetStats(false);
            histograms[pid][i]->GetXaxis()->SetLabelSize(0.04); // Increase font size for axes labels
            histograms[pid][i]->GetYaxis()->SetLabelSize(0.04);
        }
    }

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Set up TTreeReaderValues for necessary branches
    TTreeReaderValue<int> track_sector_6(mcReader, "track_sector_6");
    TTreeReaderValue<double> mc_p(mcReader, "mc_p");
    TTreeReaderValue<double> p(mcReader, "p");
    TTreeReaderValue<double> mc_theta(mcReader, "mc_theta");
    TTreeReaderValue<double> theta(mcReader, "theta");
    TTreeReaderValue<double> traj_x_6(mcReader, "traj_x_6");
    TTreeReaderValue<double> traj_y_6(mcReader, "traj_y_6");
    TTreeReaderValue<double> traj_z_6(mcReader, "traj_z_6");
    TTreeReaderValue<int> pid(mcReader, "particle_pid");

    // Edge variables for FD fiducial cuts
    TTreeReaderValue<double> edge_6(mcReader, "traj_edge_6");
    TTreeReaderValue<double> edge_18(mcReader, "traj_edge_18");
    TTreeReaderValue<double> edge_36(mcReader, "traj_edge_36");

    // Loop over events
    // for (int i = 0; i < 1e7; ++i) {
    //     mcReader.Next();
    while (mcReader.Next()) {
        if (!is_fd_track(*track_sector_6)) continue;
        if (!dc_fiducial(*edge_6, *edge_18, *edge_36, *pid)) continue;
        double delta_theta = *mc_theta - *theta;
        double theta_dc_1 = calculate_theta(*traj_x_6, *traj_y_6, *traj_z_6);

        // Check if the current particle type is one of interest and if the track is below the curve
        // if (histograms.find(*pid) != histograms.end() && !is_above_theta_dc_curve(*p, theta_dc_1)) {
        if (histograms.find(*pid) != histograms.end() ) {
            for (size_t i = 0; i < theta_bins.size(); ++i) {
                if (*theta >= theta_bins[i].first && *theta < theta_bins[i].second) {
                    histograms[*pid][i]->Fill(*p, delta_theta);
                    break;
                }
            }
        }
    }

    // Save the histograms into a canvas
    for (const auto& entry : histograms) {
        int pid = entry.first;
        const std::string& particle_name = std::get<0>(particle_types[pid]);

        TCanvas* c_deltatheta = new TCanvas(("c_deltatheta_" + particle_name).c_str(), ("Delta #theta Distributions: " + dataset + ", " + particle_name).c_str(), 2000, 1200);
        c_deltatheta->Divide(6, 4);  // 20 subplots

        std::vector<TF1*> fit_deltatheta(theta_bins.size());
        std::vector<double> A_values(theta_bins.size());
        std::vector<double> A_errors(theta_bins.size());
        std::vector<double> B_values(theta_bins.size());
        std::vector<double> B_errors(theta_bins.size());
        std::vector<double> C_values(theta_bins.size());
        std::vector<double> C_errors(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            // Ensure we are drawing on the correct pad
            c_deltatheta->cd(i + 1);
            gPad->SetMargin(0.15, 0.15, 0.20, 0.1);  // Left, right, bottom, top margins
            gPad->SetLogz();

            // Create profile histograms
            TProfile* prof_deltatheta = histograms[pid][i]->ProfileX();

            // Set the range of the profile to start from 0.5
            prof_deltatheta->GetXaxis()->SetRangeUser(0.5, std::get<2>(particle_types[pid]));

            // Fit the profiles with appropriate functions
            fit_deltatheta[i] = new TF1(("fit_deltatheta_" + std::to_string(i)).c_str(), "[0] + [1]/x + [2]/x^2", 0.5, std::get<2>(particle_types[pid]));
            prof_deltatheta->Fit(fit_deltatheta[i], "Q"); // Silent fit

            // Set the range of the fit function for plotting
            fit_deltatheta[i]->SetRange(0.5, std::get<2>(particle_types[pid]));

            // Store the fit parameters
            A_values[i] = fit_deltatheta[i]->GetParameter(0);
            A_errors[i] = fit_deltatheta[i]->GetParError(0);
            B_values[i] = fit_deltatheta[i]->GetParameter(1);
            B_errors[i] = fit_deltatheta[i]->GetParError(1);
            C_values[i] = fit_deltatheta[i]->GetParameter(2);
            C_errors[i] = fit_deltatheta[i]->GetParError(2);

            histograms[pid][i]->Draw("COLZ");
            prof_deltatheta->Draw("same");  // Draw the profile to show the fit line
            fit_deltatheta[i]->Draw("same");  // Draw the fit on top of the profile
        }

        // Add centered text "dataset, FD" on the canvas
        c_deltatheta->cd();  // Switch to the main canvas (not any specific pad)
        TLatex latex;
        latex.SetNDC();  // Use normalized coordinates (0,0) to (1,1)
        latex.SetTextSize(0.035);  // Set the text size
        latex.DrawLatex(0.425, 0.5, (dataset + ", FD").c_str());  // Add text in the center

        // Save the canvas
        c_deltatheta->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/delta_theta_distributions_" + particle_name + ".png").c_str());

        // Use the new modular function for the fitted parameters
        plot_and_fit_parameters(theta_bins, A_values, A_errors, B_values, B_errors, C_values, C_errors, particle_name, dataset, "#theta");

        // Clean up memory
        for (size_t i = 0; i < theta_bins.size(); ++i) {
            delete fit_deltatheta[i];
            delete histograms[pid][i];
        }
        delete c_deltatheta;
    }
}

void energy_loss_distributions_delta_phi_fd(TTreeReader& mcReader, const std::string& dataset) {
    // Particle types and their corresponding LaTeX names and x-axis ranges
    std::map<int, std::tuple<std::string, double, double>> particle_types = {
        {2212, {"p", 0.0, 5.0}}
    };

    // // Define theta bins
    std::vector<std::pair<double, double>> theta_bins;
    if (dataset == "rga_fa18_out") {
        // Outbending case: 24 equally spaced bins between 8 and 44 degrees
        theta_bins = {
            {8.0, 9.5}, {9.5, 11.0}, {11.0, 12.5}, {12.5, 14.0},
            {14.0, 15.5}, {15.5, 17.0}, {17.0, 18.5}, {18.5, 20.0},
            {20.0, 21.5}, {21.5, 23.0}, {23.0, 24.5}, {24.5, 26.0},
            {26.0, 27.5}, {27.5, 29.0}, {29.0, 30.5}, {30.5, 32.0},
            {32.0, 33.5}, {33.5, 35.0}, {35.0, 36.5}, {36.5, 38.0},
            {38.0, 39.5}, {39.5, 41.0}, {41.0, 42.5}, {42.5, 44.0}
        };
    } else {
        // Inbending case (default): existing bins
        theta_bins = {
            {5.0, 6.3043}, {6.3043, 7.6087}, {7.6087, 8.9130},
            {8.9130, 10.2174}, {10.2174, 11.5217}, {11.5217, 12.8261},
            {12.8261, 14.1304}, {14.1304, 15.4348}, {15.4348, 16.7391},
            {16.7391, 18.0435}, {18.0435, 19.3478}, {19.3478, 20.6522},
            {20.6522, 21.9565}, {21.9565, 23.2609}, {23.2609, 24.5652},
            {24.5652, 25.8696}, {25.8696, 27.1739}, {27.1739, 28.4783},
            {28.4783, 29.7826}, {29.7826, 31.0870}, {31.0870, 32.3913},
            {32.3913, 33.6957}, {33.6957, 35.0}, {35.0, 42.0}
        };
    }

    // Create histograms for each particle type and theta bin
    std::map<int, std::vector<TH2D*>> histograms;
    for (const auto& particle : particle_types) {
        int pid = particle.first;
        const std::string& particle_name = std::get<0>(particle.second);
        double xMin = std::get<1>(particle.second);
        double xMax = std::get<2>(particle.second);

        histograms[pid].resize(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            std::string bin_label = TString::Format("#theta [%.1f, %.1f]", theta_bins[i].first, theta_bins[i].second).Data();

            histograms[pid][i] = new TH2D(
                ("h_deltaphi_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label.c_str(),
                75, xMin, xMax, 75, -1, 1);

            // Set axis labels
            histograms[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms[pid][i]->GetYaxis()->SetTitle("#Deltaphi");

            histograms[pid][i]->SetStats(false);
            histograms[pid][i]->GetXaxis()->SetLabelSize(0.04); // Increase font size for axes labels
            histograms[pid][i]->GetYaxis()->SetLabelSize(0.04);
        }
    }

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Set up TTreeReaderValues for necessary branches
    TTreeReaderValue<int> track_sector_6(mcReader, "track_sector_6");
    TTreeReaderValue<double> mc_p(mcReader, "mc_p");
    TTreeReaderValue<double> p(mcReader, "p");
    TTreeReaderValue<double> mc_theta(mcReader, "mc_theta");
    TTreeReaderValue<double> theta(mcReader, "theta");
    TTreeReaderValue<double> mc_phi(mcReader, "mc_phi");
    TTreeReaderValue<double> phi(mcReader, "phi");
    TTreeReaderValue<double> traj_x_6(mcReader, "traj_x_6");
    TTreeReaderValue<double> traj_y_6(mcReader, "traj_y_6");
    TTreeReaderValue<double> traj_z_6(mcReader, "traj_z_6");
    TTreeReaderValue<int> pid(mcReader, "particle_pid");

    // Edge variables for FD fiducial cuts
    TTreeReaderValue<double> edge_6(mcReader, "traj_edge_6");
    TTreeReaderValue<double> edge_18(mcReader, "traj_edge_18");
    TTreeReaderValue<double> edge_36(mcReader, "traj_edge_36");

    // Loop over events
    // for (int i = 0; i < 1e7; ++i) {
    //     mcReader.Next();
    while (mcReader.Next()) {
        if (!is_fd_track(*track_sector_6)) continue;
        if (!dc_fiducial(*edge_6, *edge_18, *edge_36, *pid)) continue;
        double delta_phi = *mc_phi - *phi;
        double theta_dc_1 = calculate_theta(*traj_x_6, *traj_y_6, *traj_z_6);

        // Check if the current particle type is one of interest and if the track is below the curve
        // if (histograms.find(*pid) != histograms.end() && !is_above_theta_dc_curve(*p, theta_dc_1)) {
        if (histograms.find(*pid) != histograms.end() ) {
            for (size_t i = 0; i < theta_bins.size(); ++i) {
                if (*theta >= theta_bins[i].first && *theta < theta_bins[i].second) {
                    histograms[*pid][i]->Fill(*p, delta_phi);
                    break;
                }
            }
        }
    }

    // Save the histograms into a canvas
    for (const auto& entry : histograms) {
        int pid = entry.first;
        const std::string& particle_name = std::get<0>(particle_types[pid]);

        TCanvas* c_deltaphi = new TCanvas(("c_deltaphi_" + particle_name).c_str(), ("Delta #phi Distributions: " + dataset + ", " + particle_name).c_str(), 2000, 1200);
        c_deltaphi->Divide(6, 4);  // 20 subplots

        std::vector<TF1*> fit_deltaphi(theta_bins.size());
        std::vector<double> A_values(theta_bins.size());
        std::vector<double> A_errors(theta_bins.size());
        std::vector<double> B_values(theta_bins.size());
        std::vector<double> B_errors(theta_bins.size());
        std::vector<double> C_values(theta_bins.size());
        std::vector<double> C_errors(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            // Ensure we are drawing on the correct pad
            c_deltaphi->cd(i + 1);
            gPad->SetMargin(0.15, 0.15, 0.20, 0.1);  // Left, right, bottom, top margins
            gPad->SetLogz();

            // Create profile histograms
            TProfile* prof_deltaphi = histograms[pid][i]->ProfileX();

            // Set the range of the profile to start from 0.5
            prof_deltaphi->GetXaxis()->SetRangeUser(0.5, std::get<2>(particle_types[pid]));

            // Fit the profiles with appropriate functions
            fit_deltaphi[i] = new TF1(("fit_deltaphi_" + std::to_string(i)).c_str(), "[0] + [1]/x + [2]/x^2", 0.5, std::get<2>(particle_types[pid]));
            prof_deltaphi->Fit(fit_deltaphi[i], "Q"); // Silent fit

            // Set the range of the fit function for plotting
            fit_deltaphi[i]->SetRange(0.5, std::get<2>(particle_types[pid]));

            // Store the fit parameters
            A_values[i] = fit_deltaphi[i]->GetParameter(0);
            A_errors[i] = fit_deltaphi[i]->GetParError(0);
            B_values[i] = fit_deltaphi[i]->GetParameter(1);
            B_errors[i] = fit_deltaphi[i]->GetParError(1);
            C_values[i] = fit_deltaphi[i]->GetParameter(2);
            C_errors[i] = fit_deltaphi[i]->GetParError(2);

            histograms[pid][i]->Draw("COLZ");
            prof_deltaphi->Draw("same");  // Draw the profile to show the fit line
            fit_deltaphi[i]->Draw("same");  // Draw the fit on top of the profile
        }

        // Add centered text "dataset, FD" on the canvas
        c_deltaphi->cd();  // Switch to the main canvas (not any specific pad)
        TLatex latex;
        latex.SetNDC();  // Use normalized coordinates (0,0) to (1,1)
        latex.SetTextSize(0.035);  // Set the text size
        latex.DrawLatex(0.425, 0.5, (dataset + ", FD").c_str());  // Add text in the center

        // Save the canvas
        c_deltaphi->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/delta_phi_distributions_" + particle_name + ".png").c_str());

        // Use the new modular function for the fitted parameters
        plot_and_fit_parameters(theta_bins, A_values, A_errors, B_values, B_errors, C_values, C_errors, particle_name, dataset, "#phi");

        // Clean up memory
        for (size_t i = 0; i < theta_bins.size(); ++i) {
            delete fit_deltaphi[i];
            delete histograms[pid][i];
        }
        delete c_deltaphi;
    }
}

void plot_and_fit_parameters_cd(const std::vector<std::pair<double, double>>& theta_bins,
                             const std::vector<double>& A_values,
                             const std::vector<double>& A_errors,
                             const std::vector<double>& B_values,
                             const std::vector<double>& B_errors,
                             const std::vector<double>& C_values,
                             const std::vector<double>& C_errors,
                             const std::string& particle_name,
                             const std::string& dataset,
                             const std::string& prefix) {
    // Create a new canvas for the fitted parameters
    TCanvas* c_fit_params = new TCanvas(("c_fit_params_" + particle_name).c_str(), 
                                        ("Fit Parameters: " + dataset + ", " + particle_name).c_str(), 
                                        1600, 800);
    c_fit_params->Divide(3, 1);  // 1 row, 2 columns

    // Plot A(#theta)
    c_fit_params->cd(1);
    TGraphErrors* graph_A = new TGraphErrors(theta_bins.size());
    for (size_t i = 0; i < theta_bins.size(); ++i) {
        double theta_midpoint = 0.5 * (theta_bins[i].first + theta_bins[i].second);
        graph_A->SetPoint(i, theta_midpoint, A_values[i]);
        graph_A->SetPointError(i, 0.0, A_errors[i]);
    }
    
    if (prefix == "p") {
        graph_A->GetYaxis()->SetRangeUser(-0.02, 0.02);  // Set y-axis range
        graph_A->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_A->SetTitle(("A_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", CD;#theta (degrees);A_{" + prefix + "}(#theta) (GeV)").c_str());
    } else {
        graph_A->GetYaxis()->SetRangeUser(-0.5, 0.5);  // Set y-axis range
        graph_A->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_A->SetTitle(("A_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", CD;#theta (degrees);A_{" + prefix + "}(#theta)").c_str());
    }
    
    graph_A->SetMarkerStyle(20);  // Set marker style to a filled circle
    gPad->SetLeftMargin(0.2);  // Increase left margin
    graph_A->Draw("AP");

    // Fit A(#theta) to a 2nd order polynomial
    TF1* fit_A = new TF1("fit_A", "[0]+[1]*x+[2]*x*x", theta_bins.front().first, theta_bins.back().second);
    graph_A->Fit(fit_A, "Q");  // Silent fit
    fit_A->Draw("same");

    // Add fit results and chi2/ndf to the plot
    TPaveText* pt_A = new TPaveText(0.7, 0.75, 0.9, 0.9, "NDC");
    pt_A->AddText(Form("p0 = %.7f", fit_A->GetParameter(0)));
    pt_A->AddText(Form("p1 = %.7f", fit_A->GetParameter(1)));
    pt_A->AddText(Form("p2 = %.7f", fit_A->GetParameter(2)));
    pt_A->AddText(Form("#chi^{2}/ndf = %.3f", fit_A->GetChisquare() / fit_A->GetNDF()));
    pt_A->SetBorderSize(1);
    pt_A->SetFillColor(0);
    pt_A->Draw();

    // Plot B(#theta)
    c_fit_params->cd(2);
    TGraphErrors* graph_B = new TGraphErrors(theta_bins.size());
    for (size_t i = 0; i < theta_bins.size(); ++i) {
        double theta_midpoint = 0.5 * (theta_bins[i].first + theta_bins[i].second);
        graph_B->SetPoint(i, theta_midpoint, B_values[i]);
        graph_B->SetPointError(i, 0.0, B_errors[i]);
    }
    
    if (prefix == "p") {
        graph_B->GetYaxis()->SetRangeUser(-0.02, 0.02);  // Set y-axis range
        graph_B->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_B->SetTitle(("B_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", CD;#theta (degrees);B_{" + prefix + "}(#theta) (GeV^{2})").c_str());
    } else {
        graph_B->GetYaxis()->SetRangeUser(-0.5, 0.5);  // Set y-axis range
        graph_B->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_B->SetTitle(("B_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", CD;#theta (degrees);B_{" + prefix + "}(#theta)").c_str());
    }
    graph_B->SetMarkerStyle(20);  // Set marker style to a filled circle
    gPad->SetLeftMargin(0.2);  // Increase left margin
    graph_B->Draw("AP");

    // Fit B(#theta) to a 2nd order polynomial
    TF1* fit_B = new TF1("fit_B", "[0]+[1]*x+[2]*x*x", theta_bins.front().first, theta_bins.back().second);
    graph_B->Fit(fit_B, "Q");  // Silent fit
    fit_B->Draw("same");

    // Add fit results and chi2/ndf to the plot
    TPaveText* pt_B = new TPaveText(0.7, 0.75, 0.9, 0.9, "NDC");
    pt_B->AddText(Form("p0 = %.7f", fit_B->GetParameter(0)));
    pt_B->AddText(Form("p1 = %.7f", fit_B->GetParameter(1)));
    pt_B->AddText(Form("p2 = %.7f", fit_B->GetParameter(2)));
    pt_B->AddText(Form("#chi^{2}/ndf = %.3f", fit_B->GetChisquare() / fit_B->GetNDF()));
    pt_B->SetBorderSize(1);
    pt_B->SetFillColor(0);
    pt_B->Draw();

    // Plot C(#theta)
    c_fit_params->cd(3);
    TGraphErrors* graph_C = new TGraphErrors(theta_bins.size());
    for (size_t i = 0; i < theta_bins.size(); ++i) {
        double theta_midpoint = 0.5 * (theta_bins[i].first + theta_bins[i].second);
        graph_C->SetPoint(i, theta_midpoint, C_values[i]);
        graph_C->SetPointError(i, 0.0, C_errors[i]);
    }
    if (prefix == "p") {
        graph_C->GetYaxis()->SetRangeUser(-0.02, 0.02);  // Set y-axis range
        graph_C->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_C->SetTitle(("C_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", CD;#theta (degrees);C_{" + prefix + "}(#theta) (GeV^{3})").c_str());
    } else {
        graph_C->GetYaxis()->SetRangeUser(-0.5, 0.5);  // Set y-axis range
        graph_C->GetXaxis()->SetRangeUser(5, 40);  // Set x-axis range
        graph_C->SetTitle(("C_{" + prefix + "}, #Delta" + prefix + ", " + dataset +", CD;#theta (degrees);C_{" + prefix + "}(#theta)").c_str());
    }
    graph_C->SetMarkerStyle(20);  // Set marker style to a filled circle
    gPad->SetLeftMargin(0.2);  // Increase left margin
    graph_C->Draw("AP");

    // Fit C(#theta) to a 2nd order polynomial
    TF1* fit_C = new TF1("fit_C", "[0]+[1]*x+[2]*x*x", theta_bins.front().first, theta_bins.back().second);
    graph_C->Fit(fit_C, "Q");  // Silent fit
    fit_C->Draw("same");

    // Add fit results and chi2/ndf to the plot
    TPaveText* pt_C = new TPaveText(0.7, 0.75, 0.9, 0.9, "NDC");
    pt_C->AddText(Form("p0 = %.7f", fit_C->GetParameter(0)));
    pt_C->AddText(Form("p1 = %.7f", fit_C->GetParameter(1)));
    pt_C->AddText(Form("p2 = %.7f", fit_C->GetParameter(2)));
    pt_C->AddText(Form("#chi^{2}/ndf = %.3f", fit_C->GetChisquare() / fit_C->GetNDF()));
    pt_C->SetBorderSize(1);
    pt_C->SetFillColor(0);
    pt_C->Draw();

    // Print out the functional form of A(theta) in LaTeX format
    std::cout << "A_" << prefix << "(\\theta) = ";
    for (int i = 0; i <= 2; ++i) {
        double coeff = fit_A->GetParameter(i);
        if (i == 0) {
            std::cout << Form("%.7f", coeff);
        } else if (i == 1) {
            std::cout << Form(" %+.7f\\theta", coeff);
        } else {
            std::cout << Form(" %+.7f\\theta^%d", coeff, i);
        }
    }
    std::cout << std::endl;

    // Print out the functional form of B(theta) in LaTeX format
    std::cout << "B_" << prefix << "(\\theta) = ";
    for (int i = 0; i <= 2; ++i) {
        double coeff = fit_B->GetParameter(i);
        if (i == 0) {
            std::cout << Form("%.8f", coeff);
        } else if (i == 1) {
            std::cout << Form(" %+.8f\\theta", coeff);
        } else {
            std::cout << Form(" %+.8f\\theta^%d", coeff, i);
        }
    }
    std::cout << std::endl;

    // Print out the functional form of C(theta) in LaTeX format
    std::cout << "C_" << prefix << "(\\theta) = ";
    for (int i = 0; i <= 2; ++i) {
        double coeff = fit_C->GetParameter(i);
        if (i == 0) {
            std::cout << Form("%.8f", coeff);
        } else if (i == 1) {
            std::cout << Form(" %+.8f\\theta", coeff);
        } else {
            std::cout << Form(" %+.8f\\theta^%d", coeff, i);
        }
    }
    std::cout << std::endl;

    // Save the fit parameters canvas
    c_fit_params->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/cd_fit_params_" + prefix + "_" + particle_name + ".png").c_str());

    // Clean up memory
    delete graph_A;
    delete graph_B;
    delete graph_C;
    delete fit_A;
    delete fit_B;
    delete fit_C;
    delete pt_A;
    delete pt_B;
    delete pt_C;
    delete c_fit_params;
}

void energy_loss_distributions_delta_p_cd(TTreeReader& mcReader, const std::string& dataset) {
    // Particle types and their corresponding LaTeX names and x-axis ranges
    std::map<int, std::tuple<std::string, double, double>> particle_types = {
        {2212, {"p", 0.0, 3}}
    };

    // // Define theta bins
    std::vector<std::pair<double, double>> theta_bins = {
        {33.0, 34.542}, {34.542, 36.083}, {36.083, 37.625}, {37.625, 39.167},
        {39.167, 40.708}, {40.708, 42.250}, {42.250, 43.792}, {43.792, 45.333},
        {45.333, 46.875}, {46.875, 48.417}, {48.417, 49.958}, {49.958, 51.500},
        {51.500, 53.042}, {53.042, 54.583}, {54.583, 56.125}, {56.125, 57.667},
        {57.667, 59.208}, {59.208, 60.750}, {60.750, 62.292}, {62.292, 63.833},
        {63.833, 65.375}, {65.375, 66.917}, {66.917, 68.458}, {68.458, 70.0}
    };

    // Create histograms for each particle type and theta bin
    std::map<int, std::vector<TH2D*>> histograms;
    for (const auto& particle : particle_types) {
        int pid = particle.first;
        const std::string& particle_name = std::get<0>(particle.second);
        double xMin = std::get<1>(particle.second);
        double xMax = std::get<2>(particle.second);

        histograms[pid].resize(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            std::string bin_label = TString::Format("#theta [%.1f, %.1f]", theta_bins[i].first, theta_bins[i].second).Data();

            histograms[pid][i] = new TH2D(
                ("h_deltap_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label.c_str(),
                50, xMin, xMax, 50, -0.2, 0.2);

            // Set axis labels
            histograms[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms[pid][i]->GetYaxis()->SetTitle("#Deltap");

            histograms[pid][i]->SetStats(false);
            histograms[pid][i]->GetXaxis()->SetLabelSize(0.04); // Increase font size for axes labels
            histograms[pid][i]->GetYaxis()->SetLabelSize(0.04);
        }
    }

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Set up TTreeReaderValues for necessary branches
    TTreeReaderValue<int> track_sector_5(mcReader, "track_sector_5");
    TTreeReaderValue<double> mc_p(mcReader, "mc_p");
    TTreeReaderValue<double> p(mcReader, "p");
    TTreeReaderValue<double> mc_theta(mcReader, "mc_theta");
    TTreeReaderValue<double> theta(mcReader, "theta");
    TTreeReaderValue<double> mc_phi(mcReader, "mc_phi");
    TTreeReaderValue<double> phi(mcReader, "phi");
    TTreeReaderValue<int> pid(mcReader, "particle_pid");

    // Edge variables for FD fiducial cuts
    TTreeReaderValue<double> traj_edge_1(mcReader, "traj_edge_1");
    TTreeReaderValue<double> traj_edge_3(mcReader, "traj_edge_3");
    TTreeReaderValue<double> traj_edge_5(mcReader, "traj_edge_5");
    TTreeReaderValue<double> traj_edge_7(mcReader, "traj_edge_7");
    TTreeReaderValue<double> traj_edge_12(mcReader, "traj_edge_12");

    // Loop over events
    // for (int i = 0; i < 1e7; ++i) {
    //     mcReader.Next();
    while (mcReader.Next()) {
        if (!is_cd_track(*track_sector_5)) continue;
        if (!cvt_fiducial(*traj_edge_1, *traj_edge_3, *traj_edge_5, *traj_edge_7, *traj_edge_12)) continue;
        double delta_p = *mc_p - *p;

        // Check if the current particle type is one of interest and if the track is below the curve
        if (histograms.find(*pid) != histograms.end() ) {
            for (size_t i = 0; i < theta_bins.size(); ++i) {
                if (*theta >= theta_bins[i].first && *theta < theta_bins[i].second) {
                    histograms[*pid][i]->Fill(*p, delta_p);
                    break;
                }
            }
        }
    }

    // Save the histograms into a canvas
    for (const auto& entry : histograms) {
        int pid = entry.first;
        const std::string& particle_name = std::get<0>(particle_types[pid]);

        TCanvas* c_deltap = new TCanvas(("c_deltap_" + particle_name).c_str(), ("Delta p Distributions: " + dataset + ", " + particle_name).c_str(), 2000, 1200);
        c_deltap->Divide(6, 4);  // 20 subplots

        std::vector<TF1*> fit_deltap(theta_bins.size());
        std::vector<double> A_values(theta_bins.size());
        std::vector<double> A_errors(theta_bins.size());
        std::vector<double> B_values(theta_bins.size());
        std::vector<double> B_errors(theta_bins.size());
        std::vector<double> C_values(theta_bins.size());
        std::vector<double> C_errors(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            // Ensure we are drawing on the correct pad
            c_deltap->cd(i + 1);
            gPad->SetMargin(0.15, 0.15, 0.20, 0.1);  // Left, right, bottom, top margins
            gPad->SetLogz();

            // Create profile histograms
            TProfile* prof_deltap = histograms[pid][i]->ProfileX();

            // Set the range of the profile to start from 0.5
            prof_deltap->GetXaxis()->SetRangeUser(0.3, std::get<2>(particle_types[pid]));

            // Fit the profiles with appropriate functions
            fit_deltap[i] = new TF1(("fit_deltap_" + std::to_string(i)).c_str(), "[0] + [1]/x + [2]/x^2", 0.3, std::get<2>(particle_types[pid]));
            prof_deltap->Fit(fit_deltap[i], "Q"); // Silent fit

            // Set the range of the fit function for plotting
            fit_deltap[i]->SetRange(0.3, std::get<2>(particle_types[pid]));

            // Store the fit parameters
            A_values[i] = fit_deltap[i]->GetParameter(0);
            A_errors[i] = fit_deltap[i]->GetParError(0);
            B_values[i] = fit_deltap[i]->GetParameter(1);
            B_errors[i] = fit_deltap[i]->GetParError(1);
            C_values[i] = fit_deltap[i]->GetParameter(2);
            C_errors[i] = fit_deltap[i]->GetParError(2);

            histograms[pid][i]->Draw("COLZ");
            prof_deltap->Draw("same");  // Draw the profile to show the fit line
            fit_deltap[i]->Draw("same");  // Draw the fit on top of the profile
        }

        // Add centered text "dataset, FD" on the canvas
        c_deltap->cd();  // Switch to the main canvas (not any specific pad)
        TLatex latex;
        latex.SetNDC();  // Use normalized coordinates (0,0) to (1,1)
        latex.SetTextSize(0.035);  // Set the text size
        latex.DrawLatex(0.425, 0.5, (dataset + ", CD").c_str());  // Add text in the center

        // Save the canvas
        c_deltap->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/cd_delta_p_distributions_" + particle_name + ".png").c_str());

        // Use the new modular function for the fitted parameters
        plot_and_fit_parameters_cd(theta_bins, A_values, A_errors, B_values, B_errors, C_values, C_errors, particle_name, dataset, "p");

        // Clean up memory
        for (size_t i = 0; i < theta_bins.size(); ++i) {
            delete fit_deltap[i];
            delete histograms[pid][i];
        }
        delete c_deltap;
    }
}

void energy_loss_distributions_delta_theta_cd(TTreeReader& mcReader, const std::string& dataset) {
    // Particle types and their corresponding LaTeX names and x-axis ranges
    std::map<int, std::tuple<std::string, double, double>> particle_types = {
        {2212, {"p", 0.0, 5.0}}
    };

    // // Define theta bins
    std::vector<std::pair<double, double>> theta_bins = {
        {20.0, 22.9167}, {22.9167, 25.8333}, {25.8333, 28.75}, {28.75, 31.6667},
        {31.6667, 34.5833}, {34.5833, 37.5}, {37.5, 40.4167}, {40.4167, 43.3333},
        {43.3333, 46.25}, {46.25, 49.1667}, {49.1667, 52.0833}, {52.0833, 55.0},
        {55.0, 57.9167}, {57.9167, 60.8333}, {60.8333, 63.75}, {63.75, 66.6667},
        {66.6667, 69.5833}, {69.5833, 72.5}, {72.5, 75.4167}, {75.4167, 78.3333},
        {78.3333, 81.25}, {81.25, 84.1667}, {84.1667, 87.0833}, {87.0833, 90.0}
    };

    // Create histograms for each particle type and theta bin
    std::map<int, std::vector<TH2D*>> histograms;
    for (const auto& particle : particle_types) {
        int pid = particle.first;
        const std::string& particle_name = std::get<0>(particle.second);
        double xMin = std::get<1>(particle.second);
        double xMax = std::get<2>(particle.second);

        histograms[pid].resize(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            std::string bin_label = TString::Format("#theta [%.1f, %.1f]", theta_bins[i].first, theta_bins[i].second).Data();

            histograms[pid][i] = new TH2D(
                ("h_deltatheta_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label.c_str(),
                75, xMin, xMax, 75, -1, 1);

            // Set axis labels
            histograms[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms[pid][i]->GetYaxis()->SetTitle("#Delta#theta");

            histograms[pid][i]->SetStats(false);
            histograms[pid][i]->GetXaxis()->SetLabelSize(0.04); // Increase font size for axes labels
            histograms[pid][i]->GetYaxis()->SetLabelSize(0.04);
        }
    }

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Set up TTreeReaderValues for necessary branches
    TTreeReaderValue<int> track_sector_5(mcReader, "track_sector_5");
    TTreeReaderValue<double> mc_p(mcReader, "mc_p");
    TTreeReaderValue<double> p(mcReader, "p");
    TTreeReaderValue<double> mc_theta(mcReader, "mc_theta");
    TTreeReaderValue<double> theta(mcReader, "theta");
    TTreeReaderValue<double> mc_phi(mcReader, "mc_phi");
    TTreeReaderValue<double> phi(mcReader, "phi");
    TTreeReaderValue<int> pid(mcReader, "particle_pid");

    // Edge variables for FD fiducial cuts
    TTreeReaderValue<double> traj_edge_1(mcReader, "traj_edge_1");
    TTreeReaderValue<double> traj_edge_3(mcReader, "traj_edge_3");
    TTreeReaderValue<double> traj_edge_5(mcReader, "traj_edge_5");
    TTreeReaderValue<double> traj_edge_7(mcReader, "traj_edge_7");
    TTreeReaderValue<double> traj_edge_12(mcReader, "traj_edge_12");

    // Loop over events
    // for (int i = 0; i < 1e7; ++i) {
    //     mcReader.Next();
    while (mcReader.Next()) {
        if (!is_cd_track(*track_sector_5)) continue;
        if (!cvt_fiducial(*traj_edge_1, *traj_edge_3, *traj_edge_5, *traj_edge_7, *traj_edge_12)) continue;
        double delta_theta = *mc_theta - *theta;

        // Check if the current particle type is one of interest and if the track is below the curve
        if (histograms.find(*pid) != histograms.end() ) {
            for (size_t i = 0; i < theta_bins.size(); ++i) {
                if (*theta >= theta_bins[i].first && *theta < theta_bins[i].second) {
                    histograms[*pid][i]->Fill(*p, delta_theta);
                    break;
                }
            }
        }
    }

    // Save the histograms into a canvas
    for (const auto& entry : histograms) {
        int pid = entry.first;
        const std::string& particle_name = std::get<0>(particle_types[pid]);

        TCanvas* c_deltatheta = new TCanvas(("c_deltatheta_" + particle_name).c_str(), ("Delta #theta Distributions: " + dataset + ", " + particle_name).c_str(), 2000, 1200);
        c_deltatheta->Divide(6, 4);  // 20 subplots

        std::vector<TF1*> fit_deltatheta(theta_bins.size());
        std::vector<double> A_values(theta_bins.size());
        std::vector<double> A_errors(theta_bins.size());
        std::vector<double> B_values(theta_bins.size());
        std::vector<double> B_errors(theta_bins.size());
        std::vector<double> C_values(theta_bins.size());
        std::vector<double> C_errors(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            // Ensure we are drawing on the correct pad
            c_deltatheta->cd(i + 1);
            gPad->SetMargin(0.15, 0.15, 0.20, 0.1);  // Left, right, bottom, top margins
            gPad->SetLogz();

            // Create profile histograms
            TProfile* prof_deltatheta = histograms[pid][i]->ProfileX();

            // Set the range of the profile to start from 0.5
            prof_deltatheta->GetXaxis()->SetRangeUser(0.3, std::get<2>(particle_types[pid]));

            // Fit the profiles with appropriate functions
            fit_deltatheta[i] = new TF1(("fit_deltatheta_" + std::to_string(i)).c_str(), "[0] + [1]/x + [2]/x^2", 0.3, std::get<2>(particle_types[pid]));
            prof_deltatheta->Fit(fit_deltatheta[i], "Q"); // Silent fit

            // Set the range of the fit function for plotting
            fit_deltatheta[i]->SetRange(0.3, std::get<2>(particle_types[pid]));

            // Store the fit parameters
            A_values[i] = fit_deltatheta[i]->GetParameter(0);
            A_errors[i] = fit_deltatheta[i]->GetParError(0);
            B_values[i] = fit_deltatheta[i]->GetParameter(1);
            B_errors[i] = fit_deltatheta[i]->GetParError(1);
            C_values[i] = fit_deltatheta[i]->GetParameter(2);
            C_errors[i] = fit_deltatheta[i]->GetParError(2);

            histograms[pid][i]->Draw("COLZ");
            prof_deltatheta->Draw("same");  // Draw the profile to show the fit line
            fit_deltatheta[i]->Draw("same");  // Draw the fit on top of the profile
        }

        // Add centered text "dataset, FD" on the canvas
        c_deltatheta->cd();  // Switch to the main canvas (not any specific pad)
        TLatex latex;
        latex.SetNDC();  // Use normalized coordinates (0,0) to (1,1)
        latex.SetTextSize(0.035);  // Set the text size
        latex.DrawLatex(0.425, 0.5, (dataset + ", CD").c_str());  // Add text in the center

        // Save the canvas
        c_deltatheta->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/cd_delta_theta_distributions_" + particle_name + ".png").c_str());

        // Use the new modular function for the fitted parameters
        plot_and_fit_parameters_cd(theta_bins, A_values, A_errors, B_values, B_errors, C_values, C_errors, particle_name, dataset, "#theta");

        // Clean up memory
        for (size_t i = 0; i < theta_bins.size(); ++i) {
            delete fit_deltatheta[i];
            delete histograms[pid][i];
        }
        delete c_deltatheta;
    }
}

void energy_loss_distributions_delta_phi_cd(TTreeReader& mcReader, const std::string& dataset) {
    // Particle types and their corresponding LaTeX names and x-axis ranges
    std::map<int, std::tuple<std::string, double, double>> particle_types = {
        {2212, {"p", 0.0, 5.0}}
    };

    // // Define theta bins
    std::vector<std::pair<double, double>> theta_bins = {
        {20.0, 22.9167}, {22.9167, 25.8333}, {25.8333, 28.75}, {28.75, 31.6667},
        {31.6667, 34.5833}, {34.5833, 37.5}, {37.5, 40.4167}, {40.4167, 43.3333},
        {43.3333, 46.25}, {46.25, 49.1667}, {49.1667, 52.0833}, {52.0833, 55.0},
        {55.0, 57.9167}, {57.9167, 60.8333}, {60.8333, 63.75}, {63.75, 66.6667},
        {66.6667, 69.5833}, {69.5833, 72.5}, {72.5, 75.4167}, {75.4167, 78.3333},
        {78.3333, 81.25}, {81.25, 84.1667}, {84.1667, 87.0833}, {87.0833, 90.0}
    };

    // Create histograms for each particle type and theta bin
    std::map<int, std::vector<TH2D*>> histograms;
    for (const auto& particle : particle_types) {
        int pid = particle.first;
        const std::string& particle_name = std::get<0>(particle.second);
        double xMin = std::get<1>(particle.second);
        double xMax = std::get<2>(particle.second);

        histograms[pid].resize(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            std::string bin_label = TString::Format("#theta [%.1f, %.1f]", theta_bins[i].first, theta_bins[i].second).Data();

            histograms[pid][i] = new TH2D(
                ("h_deltaphi_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label.c_str(),
                75, xMin, xMax, 75, -1, 1);

            // Set axis labels
            histograms[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms[pid][i]->GetYaxis()->SetTitle("#Deltaphi");

            histograms[pid][i]->SetStats(false);
            histograms[pid][i]->GetXaxis()->SetLabelSize(0.04); // Increase font size for axes labels
            histograms[pid][i]->GetYaxis()->SetLabelSize(0.04);
        }
    }

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Set up TTreeReaderValues for necessary branches
    TTreeReaderValue<int> track_sector_5(mcReader, "track_sector_5");
    TTreeReaderValue<double> mc_p(mcReader, "mc_p");
    TTreeReaderValue<double> p(mcReader, "p");
    TTreeReaderValue<double> mc_theta(mcReader, "mc_theta");
    TTreeReaderValue<double> theta(mcReader, "theta");
    TTreeReaderValue<double> mc_phi(mcReader, "mc_phi");
    TTreeReaderValue<double> phi(mcReader, "phi");
    TTreeReaderValue<int> pid(mcReader, "particle_pid");

    // Edge variables for FD fiducial cuts
    TTreeReaderValue<double> traj_edge_1(mcReader, "traj_edge_1");
    TTreeReaderValue<double> traj_edge_3(mcReader, "traj_edge_3");
    TTreeReaderValue<double> traj_edge_5(mcReader, "traj_edge_5");
    TTreeReaderValue<double> traj_edge_7(mcReader, "traj_edge_7");
    TTreeReaderValue<double> traj_edge_12(mcReader, "traj_edge_12");

    // Loop over events
    // for (int i = 0; i < 1e7; ++i) {
    //     mcReader.Next();
    while (mcReader.Next()) {
        if (!is_cd_track(*track_sector_5)) continue;
        if (!cvt_fiducial(*traj_edge_1, *traj_edge_3, *traj_edge_5, *traj_edge_7, *traj_edge_12)) continue;
        double delta_phi = *mc_phi - *phi;

        // Check if the current particle type is one of interest and if the track is below the curve
        if (histograms.find(*pid) != histograms.end() ) {
            for (size_t i = 0; i < theta_bins.size(); ++i) {
                if (*theta >= theta_bins[i].first && *theta < theta_bins[i].second) {
                    histograms[*pid][i]->Fill(*p, delta_phi);
                    break;
                }
            }
        }
    }

    // Save the histograms into a canvas
    for (const auto& entry : histograms) {
        int pid = entry.first;
        const std::string& particle_name = std::get<0>(particle_types[pid]);

        TCanvas* c_deltaphi = new TCanvas(("c_deltaphi_" + particle_name).c_str(), ("Delta #phi Distributions: " + dataset + ", " + particle_name).c_str(), 2000, 1200);
        c_deltaphi->Divide(6, 4);  // 20 subplots

        std::vector<TF1*> fit_deltaphi(theta_bins.size());
        std::vector<double> A_values(theta_bins.size());
        std::vector<double> A_errors(theta_bins.size());
        std::vector<double> B_values(theta_bins.size());
        std::vector<double> B_errors(theta_bins.size());
        std::vector<double> C_values(theta_bins.size());
        std::vector<double> C_errors(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            // Ensure we are drawing on the correct pad
            c_deltaphi->cd(i + 1);
            gPad->SetMargin(0.15, 0.15, 0.20, 0.1);  // Left, right, bottom, top margins
            gPad->SetLogz();

            // Create profile histograms
            TProfile* prof_deltaphi = histograms[pid][i]->ProfileX();

            // Set the range of the profile to start from 0.5
            prof_deltaphi->GetXaxis()->SetRangeUser(0.3, std::get<2>(particle_types[pid]));

            // Fit the profiles with appropriate functions
            fit_deltaphi[i] = new TF1(("fit_deltaphi_" + std::to_string(i)).c_str(), "[0] + [1]/x + [2]/x^2", 0.3, std::get<2>(particle_types[pid]));
            prof_deltaphi->Fit(fit_deltaphi[i], "Q"); // Silent fit

            // Set the range of the fit function for plotting
            fit_deltaphi[i]->SetRange(0.3, std::get<2>(particle_types[pid]));

            // Store the fit parameters
            A_values[i] = fit_deltaphi[i]->GetParameter(0);
            A_errors[i] = fit_deltaphi[i]->GetParError(0);
            B_values[i] = fit_deltaphi[i]->GetParameter(1);
            B_errors[i] = fit_deltaphi[i]->GetParError(1);
            C_values[i] = fit_deltaphi[i]->GetParameter(2);
            C_errors[i] = fit_deltaphi[i]->GetParError(2);

            histograms[pid][i]->Draw("COLZ");
            prof_deltaphi->Draw("same");  // Draw the profile to show the fit line
            fit_deltaphi[i]->Draw("same");  // Draw the fit on top of the profile
        }

        // Add centered text "dataset, FD" on the canvas
        c_deltaphi->cd();  // Switch to the main canvas (not any specific pad)
        TLatex latex;
        latex.SetNDC();  // Use normalized coordinates (0,0) to (1,1)
        latex.SetTextSize(0.035);  // Set the text size
        latex.DrawLatex(0.425, 0.5, (dataset + ", CD").c_str());  // Add text in the center

        // Save the canvas
        c_deltaphi->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/cd_delta_phi_distributions_" + particle_name + ".png").c_str());

        // Use the new modular function for the fitted parameters
        plot_and_fit_parameters_cd(theta_bins, A_values, A_errors, B_values, B_errors, C_values, C_errors, particle_name, dataset, "#phi");

        // Clean up memory
        for (size_t i = 0; i < theta_bins.size(); ++i) {
            delete fit_deltaphi[i];
            delete histograms[pid][i];
        }
        delete c_deltaphi;
    }
}

void apply_energy_loss_correction(double& p, double& theta, double& phi, const std::string& dataset, const std::string& region) {
    // Define coefficients for the dataset "rga_fa18_inb FD"
    double A_p, B_p, C_p;
    double A_theta, B_theta, C_theta;
    double A_phi, B_phi, C_phi;

    if (dataset == "rga_fa18_inb" && region == "FD") {
        // A_p, B_p, C_p
        A_p = 0.0107818 - 0.0003279 * theta + 0.0000004 * theta * theta;
        B_p = -0.01610098 + 0.00061622 * theta + 0.00000568 * theta * theta;
        C_p = 0.01298953 - 0.00065889 * theta + 0.00001058 * theta * theta;

        // A_theta, B_theta, C_theta
        A_theta = 0.0815442 - 0.0091281 * theta + 0.0001723 * theta * theta;
        B_theta = -0.19886670 + 0.02875184 * theta - 0.00066575 * theta * theta;
        C_theta = 0.14384397 - 0.01864127 * theta + 0.00041362 * theta * theta;

        // A_phi, B_phi, C_phi
        A_phi = 0.0055658 - 0.0027920 * theta - 0.0000228 * theta * theta;
        B_phi = 0.37455807 - 0.00927964 * theta + 0.00030440 * theta * theta;
        C_phi = -0.32529570 + 0.01507693 * theta - 0.00029474 * theta * theta;
    }

    // Apply corrections
    p += A_p + B_p / p + C_p / (p * p);
    theta += A_theta + B_theta / theta + C_theta / (theta * theta);
    phi += A_phi + B_phi / phi + C_phi / (phi * phi);
}

void plot_energy_loss_corrections_fd(TTreeReader& mcReader, const std::string& dataset) {
    // Define particle types and their corresponding LaTeX names and x-axis ranges
    std::map<int, std::tuple<std::string, double, double>> particle_types = {
        {2212, {"p", 0.0, 5.0}}  // Proton as an example
    };

    // Define six theta bins between 5 and 42 degrees
    std::vector<std::pair<double, double>> theta_bins = {
        {5.0, 10.0}, {10.0, 15.0}, {15.0, 20.0}, {20.0, 25.0}, {25.0, 30.0}, 
        {30.0, 40.0}
    };

    // Define histograms before and after corrections
    std::map<int, std::vector<TH2D*>> histograms_before_p;
    std::map<int, std::vector<TH2D*>> histograms_after_p;
    std::map<int, std::vector<TH2D*>> histograms_before_theta;
    std::map<int, std::vector<TH2D*>> histograms_after_theta;
    std::map<int, std::vector<TH2D*>> histograms_before_phi;
    std::map<int, std::vector<TH2D*>> histograms_after_phi;

    for (const auto& [pid, particle_info] : particle_types) {
        const auto& [particle_name, xMin, xMax] = particle_info;
        histograms_before_p[pid].resize(theta_bins.size());
        histograms_after_p[pid].resize(theta_bins.size());
        histograms_before_theta[pid].resize(theta_bins.size());
        histograms_after_theta[pid].resize(theta_bins.size());
        histograms_before_phi[pid].resize(theta_bins.size());
        histograms_after_phi[pid].resize(theta_bins.size());

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            std::string bin_label_before_p = TString::Format("#theta [%.1f, %.1f] Before corrections", theta_bins[i].first, theta_bins[i].second).Data();
            std::string bin_label_after_p = TString::Format("#theta [%.1f, %.1f] After corrections", theta_bins[i].first, theta_bins[i].second).Data();
            std::string bin_label_before_theta = TString::Format("#theta [%.1f, %.1f] Before corrections", theta_bins[i].first, theta_bins[i].second).Data();
            std::string bin_label_after_theta = TString::Format("#theta [%.1f, %.1f] After corrections", theta_bins[i].first, theta_bins[i].second).Data();
            std::string bin_label_before_phi = TString::Format("#theta [%.1f, %.1f] Before corrections", theta_bins[i].first, theta_bins[i].second).Data();
            std::string bin_label_after_phi = TString::Format("#theta [%.1f, %.1f] After corrections", theta_bins[i].first, theta_bins[i].second).Data();

            histograms_before_p[pid][i] = new TH2D(
                ("h_p_before_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label_before_p.c_str(),
                75, xMin, xMax, 75, -0.05, 0.05
            );
            histograms_before_p[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms_before_p[pid][i]->GetYaxis()->SetTitle("#Deltap");
            histograms_before_p[pid][i]->GetXaxis()->SetTitleSize(0.05);
            histograms_before_p[pid][i]->GetYaxis()->SetTitleSize(0.05);
            histograms_before_p[pid][i]->GetXaxis()->SetLabelSize(0.045);
            histograms_before_p[pid][i]->GetYaxis()->SetLabelSize(0.045);

            histograms_after_p[pid][i] = new TH2D(
                ("h_p_after_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label_after_p.c_str(),
                75, xMin, xMax, 75, -0.05, 0.05
            );
            histograms_after_p[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms_after_p[pid][i]->GetYaxis()->SetTitle("#Deltap");
            histograms_after_p[pid][i]->GetXaxis()->SetTitleSize(0.05);
            histograms_after_p[pid][i]->GetYaxis()->SetTitleSize(0.05);
            histograms_after_p[pid][i]->GetXaxis()->SetLabelSize(0.045);
            histograms_after_p[pid][i]->GetYaxis()->SetLabelSize(0.045);

            histograms_before_theta[pid][i] = new TH2D(
                ("h_theta_before_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label_before_theta.c_str(),
                75, xMin, xMax, 75, -1, 1
            );
            histograms_before_theta[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms_before_theta[pid][i]->GetYaxis()->SetTitle("#Delta theta");
            histograms_before_theta[pid][i]->GetXaxis()->SetTitleSize(0.05);
            histograms_before_theta[pid][i]->GetYaxis()->SetTitleSize(0.05);
            histograms_before_theta[pid][i]->GetXaxis()->SetLabelSize(0.045);
            histograms_before_theta[pid][i]->GetYaxis()->SetLabelSize(0.045);

            histograms_after_theta[pid][i] = new TH2D(
                ("h_theta_after_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label_after_theta.c_str(),
                75, xMin, xMax, 75, -1, 1
            );
            histograms_after_theta[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms_after_theta[pid][i]->GetYaxis()->SetTitle("#Delta theta");
            histograms_after_theta[pid][i]->GetXaxis()->SetTitleSize(0.05);
            histograms_after_theta[pid][i]->GetYaxis()->SetTitleSize(0.05);
            histograms_after_theta[pid][i]->GetXaxis()->SetLabelSize(0.045);
            histograms_after_theta[pid][i]->GetYaxis()->SetLabelSize(0.045);

            histograms_before_phi[pid][i] = new TH2D(
                ("h_phi_before_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label_before_phi.c_str(),
                75, xMin, xMax, 75, -1, 1
            );
            histograms_before_phi[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms_before_phi[pid][i]->GetYaxis()->SetTitle("#Delta #phi");
            histograms_before_phi[pid][i]->GetXaxis()->SetTitleSize(0.05);
            histograms_before_phi[pid][i]->GetYaxis()->SetTitleSize(0.05);
            histograms_before_phi[pid][i]->GetXaxis()->SetLabelSize(0.045);
            histograms_before_phi[pid][i]->GetYaxis()->SetLabelSize(0.045);

            histograms_after_phi[pid][i] = new TH2D(
                ("h_phi_after_" + particle_name + "_bin" + std::to_string(i)).c_str(),
                bin_label_after_phi.c_str(),
                75, xMin, xMax, 75, -1, 1
            );
            histograms_after_phi[pid][i]->GetXaxis()->SetTitle("p (GeV)");
            histograms_after_phi[pid][i]->GetYaxis()->SetTitle("#Delta #phi");
            histograms_after_phi[pid][i]->GetXaxis()->SetTitleSize(0.05);
            histograms_after_phi[pid][i]->GetYaxis()->SetTitleSize(0.05);
            histograms_after_phi[pid][i]->GetXaxis()->SetLabelSize(0.045);
            histograms_after_phi[pid][i]->GetYaxis()->SetLabelSize(0.045);
        }
    }

    // Set the color palette to rainbow
    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(0);

    // Set up TTreeReaderValues for necessary branches
    TTreeReaderValue<int> track_sector_6(mcReader, "track_sector_6");
    TTreeReaderValue<double> mc_p(mcReader, "mc_p");
    TTreeReaderValue<double> p(mcReader, "p");
    TTreeReaderValue<double> mc_theta(mcReader, "mc_theta");
    TTreeReaderValue<double> theta(mcReader, "theta");
    TTreeReaderValue<double> mc_phi(mcReader, "mc_phi");
    TTreeReaderValue<double> phi(mcReader, "phi");
    TTreeReaderValue<int> pid(mcReader, "particle_pid");
    // Edge variables for FD fiducial cuts
    TTreeReaderValue<double> edge_6(mcReader, "traj_edge_6");
    TTreeReaderValue<double> edge_18(mcReader, "traj_edge_18");
    TTreeReaderValue<double> edge_36(mcReader, "traj_edge_36");

    // Loop over events
    // for (int i = 0; i < 1e8; ++i) {
    //     mcReader.Next();
    while (mcReader.Next()) {
        // Check if the track passes the required cuts
        if (!is_fd_track(*track_sector_6) || !dc_fiducial(*edge_6, *edge_18, *edge_36, *pid)) continue;

        double p_corr = *p, theta_corr = *theta, phi_corr = *phi;

        // Process the event only if the particle is in the defined map
        if (histograms_before_p.find(*pid) != histograms_before_p.end()) {
            // Fill the "before" histograms
            for (size_t i = 0; i < theta_bins.size(); ++i) {
                if (*theta >= theta_bins[i].first && *theta < theta_bins[i].second) {
                    histograms_before_p[*pid][i]->Fill(*p, *mc_p - *p);
                    histograms_before_theta[*pid][i]->Fill(*p, *mc_theta - *theta);
                    histograms_before_phi[*pid][i]->Fill(*p, *mc_phi - *phi);
                    break;
                }
            }

            // Apply energy loss corrections
            apply_energy_loss_correction(p_corr, theta_corr, phi_corr, dataset, "FD");

            // Fill the "after" histograms
            for (size_t i = 0; i < theta_bins.size(); ++i) {
                if (theta_corr >= theta_bins[i].first && theta_corr < theta_bins[i].second) {
                    histograms_after_p[*pid][i]->Fill(p_corr, *mc_p - p_corr);
                    histograms_after_theta[*pid][i]->Fill(p_corr, *mc_theta - theta_corr);
                    histograms_after_phi[*pid][i]->Fill(p_corr, *mc_phi - phi_corr);
                    break;
                }
            }
        }
    }

    // Create and save the canvases for before and after corrections
    for (const auto& [pid, particle_info] : particle_types) {
        const auto& [particle_name, xMin, xMax] = particle_info;

        TCanvas* c_p = new TCanvas(("c_p_" + particle_name).c_str(), "p Distributions: Before and After Corrections", 2400, 1600);
        c_p->Divide(6, 2);  // 2x6 subplots
        TCanvas* c_theta = new TCanvas(("c_theta_" + particle_name).c_str(), "Theta Distributions: Before and After Corrections", 2400, 1600);
        c_theta->Divide(6, 2);  // 2x6 subplots
        TCanvas* c_phi = new TCanvas(("c_phi_" + particle_name).c_str(), "Phi Distributions: Before and After Corrections", 2400, 1600);
        c_phi->Divide(6, 2);  // 2x6 subplots

        for (size_t i = 0; i < theta_bins.size(); ++i) {
            // p distribution
            c_p->cd(i + 1);  // Top row for "before"
            gPad->SetMargin(0.18, 0.15, 0.15, 0.1);  // Extra padding for the left margin
            gPad->SetLogz();  // Set log scale for counts
            histograms_before_p[pid][i]->Draw("COLZ");

            c_p->cd(i + 7);  // Bottom row for "after"
            gPad->SetMargin(0.18, 0.15, 0.15, 0.1);  // Extra padding for the left margin
            gPad->SetLogz();  // Set log scale for counts
            histograms_after_p[pid][i]->Draw("COLZ");

            // theta distribution
            c_theta->cd(i + 1);  // Top row for "before"
            gPad->SetMargin(0.18, 0.15, 0.15, 0.1);  // Extra padding for the left margin
            gPad->SetLogz();  // Set log scale for counts
            histograms_before_theta[pid][i]->Draw("COLZ");

            c_theta->cd(i + 7);  // Bottom row for "after"
            gPad->SetMargin(0.18, 0.15, 0.15, 0.1);  // Extra padding for the left margin
            gPad->SetLogz();  // Set log scale for counts
            histograms_after_theta[pid][i]->Draw("COLZ");

            // phi distribution
            c_phi->cd(i + 1);  // Top row for "before"
            gPad->SetMargin(0.18, 0.15, 0.15, 0.1);  // Extra padding for the left margin
            gPad->SetLogz();  // Set log scale for counts
            histograms_before_phi[pid][i]->Draw("COLZ");

            c_phi->cd(i + 7);  // Bottom row for "after"
            gPad->SetMargin(0.18, 0.15, 0.15, 0.1);  // Extra padding for the left margin
            gPad->SetLogz();  // Set log scale for counts
            histograms_after_phi[pid][i]->Draw("COLZ");
        }

        c_p->cd();
        TLatex latex_p;
        latex_p.SetNDC();
        latex_p.SetTextSize(0.04);  // Increase text size
        latex_p.DrawLatex(0.425, 0.5, (dataset + ", FD").c_str());
        c_p->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/p_distributions_before_after_" + particle_name + ".png").c_str());

        c_theta->cd();
        TLatex latex_theta;
        latex_theta.SetNDC();
        latex_theta.SetTextSize(0.04);  // Increase text size
        latex_theta.DrawLatex(0.425, 0.5, (dataset + ", FD").c_str());
        c_theta->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/theta_distributions_before_after_" + particle_name + ".png").c_str());

        c_phi->cd();
        TLatex latex_phi;
        latex_phi.SetNDC();
        latex_phi.SetTextSize(0.04);  // Increase text size
        latex_phi.DrawLatex(0.425, 0.5, (dataset + ", FD").c_str());
        c_phi->SaveAs(("output/calibration/energy_loss/" + dataset + "/distributions/phi_distributions_before_after_" + particle_name + ".png").c_str());

        // Clean up memory
        for (size_t i = 0; i < theta_bins.size(); ++i) {
            delete histograms_before_p[pid][i];
            delete histograms_after_p[pid][i];
            delete histograms_before_theta[pid][i];
            delete histograms_after_theta[pid][i];
            delete histograms_before_phi[pid][i];
            delete histograms_after_phi[pid][i];
        }
        delete c_p;
        delete c_theta;
        delete c_phi;
    }
}



// Main function to call both energy loss distribution functions
void energy_loss(TTreeReader& mcReader, const std::string& dataset) {
    // energy_loss_distributions(mcReader, dataset);

    // mcReader.Restart();
    // energy_loss_fd_distributions(mcReader, dataset);

    // mcReader.Restart();
    // energy_loss_fd_distributions_theta_dc(mcReader, dataset);

    mcReader.Restart();
    energy_loss_distributions_delta_p_fd(mcReader, dataset);

    // mcReader.Restart();
    // energy_loss_distributions_delta_theta_fd(mcReader, dataset);

    // mcReader.Restart();
    // energy_loss_distributions_delta_phi_fd(mcReader, dataset);

    // mcReader.Restart();
    // plot_energy_loss_corrections_fd(mcReader, dataset);

    // mcReader.Restart();
    // energy_loss_distributions_delta_p_cd(mcReader, dataset);
}
                           
void create_directories() {
    // Array of directories to check/create
    std::vector<std::string> directories = {
        "output/calibration/",
        "output/calibration/ft/",
        "output/calibration/cal/",
        "output/calibration/cal/fiducial/",
        "output/calibration/cal/fiducial/pcal",
        "output/calibration/cal/fiducial/ecin",
        "output/calibration/cal/fiducial/ecout",
        "output/calibration/cc/",
        "output/calibration/dc/",
        "output/calibration/dc/positions/",
        "output/calibration/dc/determination",
        "output/calibration/cvt/chi2pid",
        "output/calibration/cvt/determination",
        "output/calibration/cvt/positions",
        "output/calibration/energy_loss/",
        "output/calibration/energy_loss/rga_fa18_inb/",
        "output/calibration/energy_loss/rga_fa18_out/",
        "output/calibration/energy_loss/rga_sp19_inb/",
        "output/calibration/energy_loss/rgc_su22_inb/",
        "output/calibration/energy_loss/rga_fa18_inb/distributions/",
        "output/calibration/energy_loss/rga_fa18_out/distributions/",
        "output/calibration/energy_loss/rga_sp19_inb/distributions/",
        "output/calibration/energy_loss/rgc_su22_inb/distributions/"
    };

    // Iterate through each directory and create if it doesn't exist
    for (const auto& dir : directories) {
        if (!gSystem->AccessPathName(dir.c_str())) {
            std::cout << "Directory " << dir << " already exists." << std::endl;
        } else {
            if (gSystem->mkdir(dir.c_str(), true) == 0) {
                std::cout << "Directory " << dir << " created successfully." << std::endl;
            } else {
                std::cerr << "Error creating directory " << dir << std::endl;
            }
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << 
        	" <data_file.root> [<mc_file.root>]" << std::endl;
        return 1;
    }

    // Check and create directories before calling any plotting functions
    create_directories();

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
    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();

    // plot_ft_xy_energy(dataReader, mcReader);
    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // plot_ft_hit_position(dataReader, mcReader);

    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // plot_pcal_fiducial_determination(dataReader, mcReader);
    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // plot_ecin_fiducial_determination(dataReader, mcReader);
    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // plot_ecout_fiducial_determination(dataReader, mcReader);

    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // plot_cal_hit_position(dataReader, mcReader);

    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // dc_fiducial_determination(dataReader, mcReader);

    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // plot_dc_hit_position(dataReader, mcReader);

    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // cvt_fiducial_determination(dataReader, mcReader);

    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // plot_cvt_hit_position(dataReader, mcReader);

    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // plot_chi2pid_cd(dataReader, mcReader);

    dataReader.Restart();
    if (mcReader) mcReader->Restart();
    if (mcReader) energy_loss(*mcReader, "rga_fa18_out");  

    // Close files
    dataFile.Close();
    if (mcFile) {
        mcFile->Close();
        delete mcFile;
    }

    return 0;
}