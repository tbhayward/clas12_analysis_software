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

bool pcal_fiducial(double lv, double lw, double lu, int sector, int strictness) {
    // Apply strictness levels for additional cuts
    switch (strictness) {
        case 1:
            if (lw < 9 || lv < 9 || lu < 19) {
                return false;
            }
            break;
        case 2:
            if (lw < 14 || lv < 14 || lu < 29) {
                return false;
            }
            break;
        case 3:
            if ((lw < 19 || lv < 19) || (lw > 264 || lv > 264) || (lu < 39)) {
                return false;
            }
            break;
        default:
            return false;
    }

    // Specific cuts for each sector
    if (sector == 1) {
        if ((lw > 74 && lw < 80) || (lw > 84 && lw < 90) || 
            (lw > 212 && lw < 218) || (lw > 224 && lw < 230)) {
            return false;
        }
    } else if (sector == 2) {
        if ((lv > 100 && lv < 118) || (lu > 112 && lu < 118)) {
            return false;
        }
    } else if (sector == 4) {
        if (lv > 230 && lv < 242) {
            return false;
        }
    } else if (sector == 6) {
        if ((lw > 174 && lw < 180) || (lw > 195 && lw < 201)) {
            return false;
        }
    }

    // If none of the cuts apply, the track is good
    return true;
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
    TTreeReaderValue<int> cal_sector(dataReader, "cal_sector");

    TTreeReaderValue<double>* mc_cal_lv_1 = nullptr;
    TTreeReaderValue<double>* mc_cal_lw_1 = nullptr;
    TTreeReaderValue<double>* mc_cal_lu_1 = nullptr;
    TTreeReaderValue<int>* mc_cal_sector = nullptr;

    if (mcReader) {
        mc_cal_lv_1 = new TTreeReaderValue<double>(*mcReader, "cal_lv_1");
        mc_cal_lw_1 = new TTreeReaderValue<double>(*mcReader, "cal_lw_1");
        mc_cal_lu_1 = new TTreeReaderValue<double>(*mcReader, "cal_lu_1");
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
                    if (pcal_fiducial(*cal_lv_1, *cal_lw_1, *cal_lu_1, *cal_sector, 1)) {
                        h_data_1->Fill(*cal_x, *cal_y);
                    }
                    if (pcal_fiducial(*cal_lv_1, *cal_lw_1, *cal_lu_1, *cal_sector, 2)) {
                        h_data_2->Fill(*cal_x, *cal_y);
                    }
                    if (pcal_fiducial(*cal_lv_1, *cal_lw_1, *cal_lu_1, *cal_sector, 3)) {
                        h_data_3->Fill(*cal_x, *cal_y);
                    }
                }
            }

            // Fill the MC histograms if available, applying the cuts
			if (mcReader) {
				while (mcReader->Next()) {
					if (**mc_particle_pid == pid && **mc_cal_x != -9999 && **mc_cal_y != -9999) {
						h_mc_0->Fill(**mc_cal_x, **mc_cal_y); // No cuts
					if (pcal_fiducial(**mc_cal_lv_1, **mc_cal_lw_1, **mc_cal_lu_1, **mc_cal_sector, 1)) {
						h_mc_1->Fill(**mc_cal_x, **mc_cal_y);
					}
					if (pcal_fiducial(**mc_cal_lv_1, **mc_cal_lw_1, **mc_cal_lu_1, **mc_cal_sector, 2)) {
						h_mc_2->Fill(**mc_cal_x, **mc_cal_y);
					}
					if (pcal_fiducial(**mc_cal_lv_1, **mc_cal_lw_1, **mc_cal_lu_1, **mc_cal_sector, 3)) {
						h_mc_3->Fill(**mc_cal_x, **mc_cal_y);
						}
					}
				}
			}
	        // Draw and save the data plots
	        TCanvas c_data_0(("c_data_0_" + particle_name + "_" + layer_name).c_str(), ("c_data_" + particle_name + "_" + layer_name).c_str(), 800, 600);
	        c_data_0.SetLogz();  // Set the z-axis to a logarithmic scale
	        h_data_0->Draw("COLZ");
	        c_data_0.SaveAs(("output/calibration/cal/" + particle_name + "_data_" + layer_name + "_cal_hit_position_strictness_0.png").c_str());

	        TCanvas c_data_1(("c_data_1_" + particle_name + "_" + layer_name).c_str(), ("c_data_" + particle_name + "_" + layer_name).c_str(), 800, 600);
	        c_data_1.SetLogz();  // Set the z-axis to a logarithmic scale
	        h_data_1->Draw("COLZ");
	        c_data_1.SaveAs(("output/calibration/cal/" + particle_name + "_data_" + layer_name + "_cal_hit_position_strictness_1.png").c_str());

	        TCanvas c_data_2(("c_data_2_" + particle_name + "_" + layer_name).c_str(), ("c_data_" + particle_name + "_" + layer_name).c_str(), 800, 600);
	        c_data_2.SetLogz();  // Set the z-axis to a logarithmic scale
	        h_data_2->Draw("COLZ");
	        c_data_2.SaveAs(("output/calibration/cal/" + particle_name + "_data_" + layer_name + "_cal_hit_position_strictness_2.png").c_str());

	        TCanvas c_data_3(("c_data_3_" + particle_name + "_" + layer_name).c_str(), ("c_data_" + particle_name + "_" + layer_name).c_str(), 800, 600);
	        c_data_3.SetLogz();  // Set the z-axis to a logarithmic scale
	        h_data_3->Draw("COLZ");
	        c_data_3.SaveAs(("output/calibration/cal/" + particle_name + "_data_" + layer_name + "_cal_hit_position_strictness_3.png").c_str());

	        // Draw and save the MC plots if available
	        if (mcReader) {
	            TCanvas c_mc_0(("c_mc_0_" + particle_name + "_" + layer_name).c_str(), ("c_mc_" + particle_name + "_" + layer_name).c_str(), 800, 600);
	            c_mc_0.SetLogz();  // Set the z-axis to a logarithmic scale
	            h_mc_0->Draw("COLZ");
	            c_mc_0.SaveAs(("output/calibration/cal/" + particle_name + "_mc_" + layer_name + "_cal_hit_position_strictness_0.png").c_str());

	            TCanvas c_mc_1(("c_mc_1_" + particle_name + "_" + layer_name).c_str(), ("c_mc_" + particle_name + "_" + layer_name).c_str(), 800, 600);
	            c_mc_1.SetLogz();  // Set the z-axis to a logarithmic scale
	            h_mc_1->Draw("COLZ");
	            c_mc_1.SaveAs(("output/calibration/cal/" + particle_name + "_mc_" + layer_name + "_cal_hit_position_strictness_1.png").c_str());

	            TCanvas c_mc_2(("c_mc_2_" + particle_name + "_" + layer_name).c_str(), ("c_mc_" + particle_name + "_" + layer_name).c_str(), 800, 600);
	            c_mc_2.SetLogz();  // Set the z-axis to a logarithmic scale
	            h_mc_2->Draw("COLZ");
	            c_mc_2.SaveAs(("output/calibration/cal/" + particle_name + "_mc_" + layer_name + "_cal_hit_position_strictness_2.png").c_str());

	            TCanvas c_mc_3(("c_mc_3_" + particle_name + "_" + layer_name).c_str(), ("c_mc_" + particle_name + "_" + layer_name).c_str(), 800, 600);
	            c_mc_3.SetLogz();  // Set the z-axis to a logarithmic scale
	            h_mc_3->Draw("COLZ");
	            c_mc_3.SaveAs(("output/calibration/cal/" + particle_name + "_mc_" + layer_name + "_cal_hit_position_strictness_3.png").c_str());
	        }

	        // Clean up for this layer and particle type
	        delete h_data_0;
	        delete h_data_1;
	        delete h_data_2;
	        delete h_data_3;
	        if (h_mc_0) delete h_mc_0;
	        if (h_mc_1) delete h_mc_1;
	        if (h_mc_2) delete h_mc_2;
	        if (h_mc_3) delete h_mc_3;
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
    TTreeReaderValue<int> cal_sector(dataReader, "cal_sector");

    TTreeReaderValue<double>* mc_cal_lv_1 = nullptr;
    TTreeReaderValue<double>* mc_cal_lw_1 = nullptr;
    TTreeReaderValue<double>* mc_cal_lu_1 = nullptr;
    TTreeReaderValue<int>* mc_cal_sector = nullptr;

    if (mcReader) {
        mc_cal_lv_1 = new TTreeReaderValue<double>(*mcReader, "cal_lv_1");
        mc_cal_lw_1 = new TTreeReaderValue<double>(*mcReader, "cal_lw_1");
        mc_cal_lu_1 = new TTreeReaderValue<double>(*mcReader, "cal_lu_1");
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
                    if (pcal_fiducial(*cal_lv_1, *cal_lw_1, *cal_lu_1, *cal_sector, 1)) {
                        h_data_1->Fill(*cal_x, *cal_y);
                    }
                    if (pcal_fiducial(*cal_lv_1, *cal_lw_1, *cal_lu_1, *cal_sector, 2)) {
                        h_data_2->Fill(*cal_x, *cal_y);
                    }
                    if (pcal_fiducial(*cal_lv_1, *cal_lw_1, *cal_lu_1, *cal_sector, 3)) {
                        h_data_3->Fill(*cal_x, *cal_y);
                    }
                }
            }
            // Fill the MC histograms if available, applying the cuts
            if (mcReader) {
                while (mcReader->Next()) {
                    if (**mc_particle_pid == pid && **mc_cal_x != -9999 && **mc_cal_y != -9999) {
                        h_mc_0->Fill(**mc_cal_x, **mc_cal_y); // No cuts
                        if (pcal_fiducial(**mc_cal_lv_1, **mc_cal_lw_1, **mc_cal_lu_1, **mc_cal_sector, 1)) {
                            h_mc_1->Fill(**mc_cal_x, **mc_cal_y);
                        }
                        if (pcal_fiducial(**mc_cal_lv_1, **mc_cal_lw_1, **mc_cal_lu_1, **mc_cal_sector, 2)) {
                            h_mc_2->Fill(**mc_cal_x, **mc_cal_y);
                        }
                        if (pcal_fiducial(**mc_cal_lv_1, **mc_cal_lw_1, **mc_cal_lu_1, **mc_cal_sector, 3)) {
                            h_mc_3->Fill(**mc_cal_x, **mc_cal_y);
                        }
                    }
                }
            }

            // Create a canvas to hold the 2x4 subplots
            TCanvas* c = new TCanvas(("c_" + particle_name + "_" + layer_name).c_str(), ("c_" + particle_name + "_" + layer_name).c_str(), 1600, 800);
            c->Divide(4, 2);

            // Draw the data plots on the top row
            c->cd(1); h_data_0->Draw("COLZ");
            c->cd(2); h_data_1->Draw("COLZ");
            c->cd(3); h_data_2->Draw("COLZ");
            c->cd(4); h_data_3->Draw("COLZ");

            // Draw the MC plots on the bottom row, if available
            if (mcReader) {
                c->cd(5); h_mc_0->Draw("COLZ");
                c->cd(6); h_mc_1->Draw("COLZ");
                c->cd(7); h_mc_2->Draw("COLZ");
                c->cd(8); h_mc_3->Draw("COLZ");
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
                           
void create_directories() {
    // Array of directories to check/create
    std::vector<std::string> directories = {
        "output/calibration/",
        "output/calibration/ft/",
        "output/calibration/cal/",
        "output/calibration/cal/fiducial/",
        "output/calibration/cc/"
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
    plot_cal_hit_position(dataReader, mcReader);
    dataReader.Restart();
    if (mcReader) mcReader->Restart();
    plot_cal_fiducial_determination(dataReader, mcReader);



    // Close files
    dataFile.Close();
    if (mcFile) {
        mcFile->Close();
        delete mcFile;
    }

    return 0;
}