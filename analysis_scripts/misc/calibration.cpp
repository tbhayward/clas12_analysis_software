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
    // RGA only so far (not RGC)
    // if (sector == 1) {
    //     if ((lw_1 > 69 && lw_1 < 96) || (lw_1 > 207 && lw_1 < 236)) {
    //         return false;
    //     }
    // } else if (sector == 2) {
    //     if ((lv_1 > 95 && lv_1 < 119) || (lu_1 > 108 && lu_1 < 126)) {
    //         return false;
    //     }
    // } else if (sector == 4) {
    //     if (lv_1 > 224 && lv_1 < 247) {
    //         return false;
    //     }
    // } else if (sector == 6) {
    //     if ((lw_1 > 169 && lw_1 < 198)) {
    //         return false;
    //     }
    // }

    // Specific cuts for each sector in ECin
    // RGA and RGC
    if (sector == 1) {
    	if (lv_4 > 72 && lv_4 < 94) {
    		return false;
    	}
    }

    // Specific cuts for each sector in ECout
    // RGC only so far (not RGA)
    if (sector == 2) {
        if (lw_7 > 68 && lw_7 < 84) {
            return false;
        }
    }
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
	int strictness) {
    std::cout << edge_6 << std::endl;
    return edge_6 > 100;
    // If none of the cuts apply, the track is good
    // return true;
}

void plot_dc_hit_position(TTreeReader& dataReader, TTreeReader* mcReader = nullptr) {
     // Define the number of bins for the histograms
    int nBins = 100;
    // Array of DC regions and their corresponding variable names
    std::vector<std::tuple<std::string, std::string, std::string, double, double>> regions = {
        {"traj_x_6", "traj_y_6", "region_1", -200, 200},
        {"traj_x_18", "traj_y_18", "region_2", -300, 300},
        {"traj_x_36", "traj_y_36", "region_3", -450, 450}
    };
    // Array of particle types (photons and electrons) and their corresponding PIDs
    std::vector<std::tuple<int, std::string>> particle_types = {
        {11, "electron"},
        // {-211, "pim"},
        // {211, "pip"},
        // {321, "kp"},
        // {-321, "km"},
        {2212, "proton"}
    };

    // Declare TTreeReaderValues for the DC edge and track variables
    TTreeReaderValue<double> traj_edge_6(dataReader, "traj_edge_6");
    TTreeReaderValue<double> traj_edge_18(dataReader, "traj_edge_18");
    TTreeReaderValue<double> traj_edge_36(dataReader, "traj_edge_36");

    TTreeReaderValue<int> track_sector_6(dataReader, "track_sector_6");
    TTreeReaderValue<double> track_chi2_6(dataReader, "track_chi2_6");
    TTreeReaderValue<int> track_ndf_6(dataReader, "track_ndf_6");

    TTreeReaderValue<double>* mc_traj_edge_6 = nullptr;
    TTreeReaderValue<double>* mc_traj_edge_18 = nullptr;
    TTreeReaderValue<double>* mc_traj_edge_36 = nullptr;

    if (mcReader) {
        mc_traj_edge_6 = new TTreeReaderValue<double>(*mcReader, "traj_edge_6");
        mc_traj_edge_18 = new TTreeReaderValue<double>(*mcReader, "traj_edge_18");
        mc_traj_edge_36 = new TTreeReaderValue<double>(*mcReader, "traj_edge_36");
    }

    // Loop over each particle type
    for (const auto& particle_type : particle_types) {
        int pid = std::get<0>(particle_type);
        std::string particle_name = std::get<1>(particle_type);

        // Loop over each DC region
        for (const auto& region : regions) {
            std::string x_branch = std::get<0>(region);
            std::string y_branch = std::get<1>(region);
            std::string region_name = std::get<2>(region);
            double xMin = std::get<3>(region);
            double xMax = std::get<4>(region);
            double yMin = xMin;  // Same as xMin
            double yMax = xMax;  // Same as xMax

            // Restart the TTreeReader to process the data from the beginning
            dataReader.Restart();
            if (mcReader) mcReader->Restart();

            // Declare TTreeReaderValues for data and MC for this region
            TTreeReaderValue<double> traj_x(dataReader, x_branch.c_str());
            TTreeReaderValue<double> traj_y(dataReader, y_branch.c_str());
            TTreeReaderValue<int> particle_pid(dataReader, "particle_pid");

            TTreeReaderValue<double>* mc_traj_x = nullptr;
            TTreeReaderValue<double>* mc_traj_y = nullptr;
            TTreeReaderValue<int>* mc_particle_pid = nullptr;

            if (mcReader) {
                mc_traj_x = new TTreeReaderValue<double>(*mcReader, x_branch.c_str());
                mc_traj_y = new TTreeReaderValue<double>(*mcReader, y_branch.c_str());
                mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
            }

            // Create histograms for data and MC for each strictness level
            TH2D* h_data_0 = new TH2D("h_data_0", ("data " + region_name + " hit position (" + particle_name + ", strictness_0)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
            TH2D* h_data_1 = new TH2D("h_data_1", ("data " + region_name + " hit position (" + particle_name + ", strictness_1)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
            TH2D* h_data_2 = new TH2D("h_data_2", ("data " + region_name + " hit position (" + particle_name + ", strictness_2)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
            TH2D* h_data_3 = new TH2D("h_data_3", ("data " + region_name + " hit position (" + particle_name + ", strictness_3)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);

            h_data_0->GetXaxis()->SetTitle(("x_{" + region_name + "}").c_str());
            h_data_0->GetYaxis()->SetTitle(("y_{" + region_name + "}").c_str());
            h_data_1->GetXaxis()->SetTitle(("x_{" + region_name + "}").c_str());
            h_data_1->GetYaxis()->SetTitle(("y_{" + region_name + "}").c_str());
            h_data_2->GetXaxis()->SetTitle(("x_{" + region_name + "}").c_str());
            h_data_2->GetYaxis()->SetTitle(("y_{" + region_name + "}").c_str());
            h_data_3->GetXaxis()->SetTitle(("x_{" + region_name + "}").c_str());
            h_data_3->GetYaxis()->SetTitle(("y_{" + region_name + "}").c_str());

            TH2D* h_mc_0 = nullptr;
            TH2D* h_mc_1 = nullptr;
            TH2D* h_mc_2 = nullptr;
            TH2D* h_mc_3 = nullptr;

            if (mcReader) {
                h_mc_0 = new TH2D("h_mc_0", ("mc " + region_name + " hit position (" + particle_name + ", strictness_0)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
                h_mc_1 = new TH2D("h_mc_1", ("mc " + region_name + " hit position (" + particle_name + ", strictness_1)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
                h_mc_2 = new TH2D("h_mc_2", ("mc " + region_name + " hit position (" + particle_name + ", strictness_2)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);
                h_mc_3 = new TH2D("h_mc_3", ("mc " + region_name + " hit position (" + particle_name + ", strictness_3)").c_str(), nBins, xMin, xMax, nBins, yMin, yMax);

                h_mc_0->GetXaxis()->SetTitle(("x_{" + region_name + "}").c_str());
                h_mc_0->GetYaxis()->SetTitle(("y_{" + region_name + "}").c_str());
                h_mc_1->GetXaxis()->SetTitle(("x_{" + region_name + "}").c_str());
                h_mc_1->GetYaxis()->SetTitle(("y_{" + region_name + "}").c_str());
                h_mc_2->GetXaxis()->SetTitle(("x_{" + region_name + "}").c_str());
                h_mc_2->GetYaxis()->SetTitle(("y_{" + region_name + "}").c_str());
                h_mc_3->GetXaxis()->SetTitle(("x_{" + region_name + "}").c_str());
                h_mc_3->GetYaxis()->SetTitle(("y_{" + region_name + "}").c_str());
            }

            // Fill the data histograms, applying the cuts
            while (dataReader.Next()) {
                if (*particle_pid == pid && *traj_x != -9999 && *traj_y != -9999) {
                    h_data_0->Fill(*traj_x, *traj_y); // No cuts
                    if (dc_fiducial(*traj_edge_6, *traj_edge_18, *traj_edge_36, 1)) {
                        h_data_1->Fill(*traj_x, *traj_y);
                    }
                    if (dc_fiducial(*traj_edge_6, *traj_edge_18, *traj_edge_36, 2)) {
                        h_data_2->Fill(*traj_x, *traj_y);
                    }
                    if (dc_fiducial(*traj_edge_6, *traj_edge_18, *traj_edge_36, 3)) {
                        h_data_3->Fill(*traj_x, *traj_y);
                    }
                }
            }
            // Fill the MC histograms if available, applying the cuts
            if (mcReader) {
                while (mcReader->Next()) {
                    if (**mc_particle_pid == pid && **mc_traj_x != -9999 && **mc_traj_y != -9999) {
                        h_mc_0->Fill(**mc_traj_x, **mc_traj_y); // No cuts
                        if (dc_fiducial(**mc_traj_edge_6, **mc_traj_edge_18, **mc_traj_edge_36, 1)) {
                            h_mc_1->Fill(**mc_traj_x, **mc_traj_y);
                        }
                        if (dc_fiducial(**mc_traj_edge_6, **mc_traj_edge_18, **mc_traj_edge_36, 2)) {
                            h_mc_2->Fill(**mc_traj_x, **mc_traj_y);
                        }
                        if (dc_fiducial(**mc_traj_edge_6, **mc_traj_edge_18, **mc_traj_edge_36, 3)) {
                            h_mc_3->Fill(**mc_traj_x, **mc_traj_y);
                        }
                    }
                }
            }

            // Create a canvas to hold the 2x4 subplots
            TCanvas* c = new TCanvas(("c_" + particle_name + "_" + region_name).c_str(), ("c_" + particle_name + "_" + region_name).c_str(), 1600, 800);
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
            c->SaveAs(("output/calibration/dc/positions/" + particle_name + "_" + region_name + "_dc_hit_position.png").c_str());

            // Clean up for this region and particle type
            delete h_data_0;
            delete h_data_1;
            delete h_data_2;
            delete h_data_3;
            if (h_mc_0) delete h_mc_0;
            if (h_mc_1) delete h_mc_1;
            if (h_mc_2) delete h_mc_2;
            if (h_mc_3) delete h_mc_3;
            delete c;
            if (mc_traj_x) delete mc_traj_x;
            if (mc_traj_y) delete mc_traj_y;
            if (mc_particle_pid) delete mc_particle_pid;
        }
    }

    // Clean up the dynamically allocated memory for edge variables
    if (mc_traj_edge_6) delete mc_traj_edge_6;
    if (mc_traj_edge_18) delete mc_traj_edge_18;
    if (mc_traj_edge_36) delete mc_traj_edge_36;
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
    std::vector<std::tuple<std::string, std::string, std::string, double, double, double, double>> regions = {
        {"traj_x_6", "traj_y_6", "region_1", 15, 160, -80, 80},
        {"traj_x_18", "traj_y_18", "region_2", 30, 240, -125, 125},
        {"traj_x_36", "traj_y_36", "region_3", 30, 400, -200, 200} 
    };

    std::vector<std::tuple<int, std::string>> particle_types = {
        {11, "electron"},
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

            TTreeReaderValue<double>* mc_traj_x = nullptr;
            TTreeReaderValue<double>* mc_traj_y = nullptr;
            TTreeReaderValue<int>* mc_particle_pid = nullptr;
            TTreeReaderValue<int>* mc_track_sector_6 = nullptr;
            TTreeReaderValue<double>* mc_track_chi2_6 = nullptr;
            TTreeReaderValue<int>* mc_track_ndf_6 = nullptr;

            if (mcReader) {
                mc_traj_x = new TTreeReaderValue<double>(*mcReader, x_branch.c_str());
                mc_traj_y = new TTreeReaderValue<double>(*mcReader, y_branch.c_str());
                mc_particle_pid = new TTreeReaderValue<int>(*mcReader, "particle_pid");
                mc_track_sector_6 = new TTreeReaderValue<int>(*mcReader, "track_sector_6");
                mc_track_chi2_6 = new TTreeReaderValue<double>(*mcReader, "track_chi2_6");
                mc_track_ndf_6 = new TTreeReaderValue<int>(*mcReader, "track_ndf_6");
            }

            // Fill data histograms
            while (dataReader.Next()) {
                if (*particle_pid == pid && *traj_x != -9999 && *traj_y != -9999 && *track_ndf_6 > 0) {
                    double chi2_ndf = *track_chi2_6 / *track_ndf_6;

                    auto rotated_coords = rotate_coordinates(*traj_x, *traj_y, *track_sector_6);
                    double traj_x_rot = rotated_coords.first;
                    double traj_y_rot = rotated_coords.second;

                    // Only fill histograms if the rotated x-coordinate is positive
                    if (traj_x_rot >= 0) {
                        h_data_sum_sector[*track_sector_6 - 1]->Fill(traj_x_rot, traj_y_rot, chi2_ndf);
                        h_data_count_sector[*track_sector_6 - 1]->Fill(traj_x_rot, traj_y_rot);
                    }
                }
            }

            // Fill MC histograms
            if (mcReader) {
                while (mcReader->Next()) {
                    if (**mc_particle_pid == pid && **mc_traj_x != -9999 && **mc_traj_y != -9999 && **mc_track_ndf_6 > 0) {
                        double mc_chi2_ndf = **mc_track_chi2_6 / **mc_track_ndf_6;

                        auto rotated_coords = rotate_coordinates(**mc_traj_x, **mc_traj_y, **mc_track_sector_6);
                        double mc_traj_x_rot = rotated_coords.first;
                        double mc_traj_y_rot = rotated_coords.second;

                        // Only fill histograms if the rotated x-coordinate is positive
                        if (mc_traj_x_rot >= 0) {
                            h_mc_sum_sector[**mc_track_sector_6 - 1]->Fill(mc_traj_x_rot, mc_traj_y_rot, mc_chi2_ndf);
                            h_mc_count_sector[**mc_track_sector_6 - 1]->Fill(mc_traj_x_rot, mc_traj_y_rot);
                        }
                    }
                }
            }

            // Normalize and save histograms
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

            // Cleanup
            delete c_region;
            if (c_mc_region) delete c_mc_region;

            for (auto& hist : h_data_sum_sector) delete hist;
            for (auto& hist : h_data_count_sector) delete hist;

            if (mcReader) {
                for (auto& hist : h_mc_sum_sector) delete hist;
                for (auto& hist : h_mc_count_sector) delete hist;
            }

            if (mc_traj_x) delete mc_traj_x;
            if (mc_traj_y) delete mc_traj_y;
            if (mc_particle_pid) delete mc_particle_pid;
            if (mc_track_sector_6) delete mc_track_sector_6;
            if (mc_track_chi2_6) delete mc_track_chi2_6;
            if (mc_track_ndf_6) delete mc_track_ndf_6;
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
        "output/calibration/cvt/chi2pid"
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

    dataReader.Restart();
    if (mcReader) mcReader->Restart();
    plot_dc_hit_position(dataReader, mcReader);

    // dataReader.Restart();
    // if (mcReader) mcReader->Restart();
    // plot_chi2pid_cd(dataReader, mcReader);



    // Close files
    dataFile.Close();
    if (mcFile) {
        mcFile->Close();
        delete mcFile;
    }

    return 0;
}