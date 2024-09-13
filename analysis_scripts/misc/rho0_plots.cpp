#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm> 
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLegendEntry.h>
#include <TAxis.h>
#include <TSystem.h>
#include <vector>
#include <string>
#include <map>
#include <TLatex.h>
#include <TText.h>  
#include <TF1.h>        
#include <TMarker.h>  
#include <TPaveText.h>


std::map<std::string, std::vector<std::vector<double>>> readAsymmetries(const std::string &filename) {
    std::map<std::string, std::vector<std::vector<double>>> asymmetryData;
    std::ifstream file(filename);
    std::string line;

    // Step 1: Read the file and populate asymmetryData
    while (std::getline(file, line)) {
        std::string::size_type pos = line.find("=");
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);  // Extract key before '='
            std::string dataStr = line.substr(pos + 1);  // Data after '='

            // Clean up the data string
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), '{'), dataStr.end());
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), '}'), dataStr.end());
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), ';'), dataStr.end());

            // Split dataStr into individual values (x, y, error)
            std::vector<std::vector<double>> values;
            std::stringstream ss(dataStr);
            std::string num;
            std::vector<double> tempVec;

            while (std::getline(ss, num, ',')) {
                if (!num.empty()) {
                    tempVec.push_back(std::stod(num));  // Convert to double and store
                }
                if (tempVec.size() == 3) {  // Once we have 3 numbers (x, y, error), push to values vector
                    values.push_back(tempVec);
                    tempVec.clear();
                }
            }

            asymmetryData[key] = values;  // Store key-value pair in map
        }
    }

    // Step 2: Debugging - Print all available keys
    std::cout << "Available keys in asymmetryData:" << std::endl;
    for (const auto& entry : asymmetryData) {
        std::cout << entry.first << std::endl;
    }

    // Step 3: Attempt to calculate the doubleratio fits
    std::string alusKeySuffix = "ALUsinphi";
    std::string allKeySuffix = "ALL";

    for (const auto &entry : asymmetryData) {
        // We assume that the baseKey should exclude the suffix part ("ALUsinphi", "ALL", etc.)
        if (entry.first.find(alusKeySuffix) != std::string::npos) {
            // Extract the base key before "ALUsinphi"
            std::string baseKey = entry.first.substr(0, entry.first.find("chi2Fits"));

            std::cout << "Processing Base Key: " << baseKey << std::endl;  // Debugging print statement

            // Construct the full keys for ALUsinphi and ALL
            std::string alusKey = baseKey + "chi2Fits" + alusKeySuffix;
            std::string allKey = baseKey + "chi2Fits" + allKeySuffix;

            // Check if the expected ALUsinphi and ALL keys exist
            bool alusExists = (asymmetryData.find(alusKey) != asymmetryData.end());
            bool allExists = (asymmetryData.find(allKey) != asymmetryData.end());

            std::cout << "Checking for ALUsinphi key: " << alusKey << " - " << (alusExists ? "Found" : "Not Found") << std::endl;
            std::cout << "Checking for ALL key: " << allKey << " - " << (allExists ? "Found" : "Not Found") << std::endl;

            // Ensure both ALUsinphi and ALL exist for this bin
            if (alusExists && allExists) {
                const auto &alusData = asymmetryData[alusKey];
                const auto &allData = asymmetryData[allKey];

                // Prepare the doubleratio data
                std::vector<std::vector<double>> doubleratioData;
                for (size_t i = 0; i < alusData.size(); ++i) {
                    double xValue = alusData[i][0];
                    double ratioValue = alusData[i][1] / allData[i][1];
                    double error = std::abs(ratioValue) * std::sqrt(
                        std::pow(alusData[i][2] / alusData[i][1], 2) + 
                        std::pow(allData[i][2] / allData[i][1], 2)
                    );
                    ratioValue = -ratioValue;  // Adjust the sign based on your requirements
                    doubleratioData.push_back({xValue, ratioValue, error});
                }

                // Store the doubleratio data
                asymmetryData[baseKey + "chi2Fitsdoubleratio"] = doubleratioData;
            } else {
                // Print a warning if either ALUsinphi or ALL data is missing
                std::cout << "Missing data for baseKey: " << baseKey << std::endl;
            }
        }
    }

    return asymmetryData;
}

// Helper functions to create graphs and set axis labels
TGraphErrors* createTGraphErrors(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& yErr,
    const int markerStyle,
    const double markerSize,
    const int color) {
    TGraphErrors* graph = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, yErr.data());
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetTitle("");
    return graph;
}

void setAxisLabelsAndRanges(
    TGraphErrors* graph,
    const std::string& xLabel,
    const std::string& yLabel,
    const std::pair<double, double>& xLimits,
    const std::pair<double, double>& yLimits) {
    graph->GetXaxis()->SetTitle(xLabel.c_str());
    graph->GetYaxis()->SetTitle(yLabel.c_str());
    graph->GetXaxis()->SetLimits(xLimits.first, xLimits.second);
    graph->GetXaxis()->SetRangeUser(xLimits.first, xLimits.second);
    graph->GetYaxis()->SetRangeUser(yLimits.first, yLimits.second);

    graph->GetXaxis()->SetLabelSize(0.04);
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetLabelSize(0.04);
    graph->GetYaxis()->SetTitleSize(0.045);
}


// Helper function to print a vector of doubles
void printVector(const std::vector<double>& vec, const std::string& vecName) {
    std::cout << vecName << ": [";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i != vec.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
}

void plotDependence(
    const std::map<std::string, std::vector<std::vector<double>>>& asymmetryData,  // Single map of all asymmetry data
    const std::vector<std::string>& prefixes,  // Vector of prefixes corresponding to the 8 different asymmetry datasets
    const std::string& xLabel, 
    const std::pair<double, double>& xLimits, 
    const std::string& outputFileName,
    const std::vector<std::string>& legendEntries) {

    TCanvas* c = new TCanvas("c", "Dependence Plots", 1200, 800);
    c->Divide(3, 2);

    std::vector<std::string> suffixes = {"ALUsinphi", "AULsinphi", "AULsin2phi", "ALL", "doubleratio", "ALLcosphi"};
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin(2#phi)}/F_{UU}",
        "F_{LL}/F_{UU}",
        "-F_{LU}^{sin#phi}/F_{LL}",
        "F_{LL}^{cos#phi}/F_{UU}"
    };

    // Define a color palette for the plots (one for each dataset)
    std::vector<int> colors = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kViolet};

    std::cout << "Available keys in asymmetryData:" << std::endl;
	for (const auto& pair : asymmetryData) {
	    std::cout << pair.first << std::endl;
	}

    // Loop over each subplot (six total)
    // for (size_t i = 0; i < suffixes.size(); ++i) {
    for (size_t i = 0; i < 1; ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.15);

        std::vector<TGraphErrors*> graphs;

        // Iterate over the 8 datasets (prefixes)
        for (size_t datasetIndex = 0; datasetIndex < prefixes.size(); ++datasetIndex) {
        	
            std::string key = prefixes[datasetIndex] + "chi2Fits" + suffixes[i];
            auto it = asymmetryData.find(key);

            if (it != asymmetryData.end()) {
                const auto& data = it->second;

                std::vector<double> x, y, yStatErr;
                for (const auto& entry : data) {
                    x.push_back(entry[0]);
                    y.push_back(entry[1]);
                    yStatErr.push_back(entry[2]);
                }
                
                // Print the vectors
				printVector(x, "x");
				printVector(y, "y");
				printVector(yStatErr, "yStatErr");

                TGraphErrors* graph = createTGraphErrors(x, y, yStatErr, 20, 0.8, colors[datasetIndex]);
                graphs.push_back(graph);

                // Set the axis labels and ranges for the first dataset
                if (datasetIndex == 0) {
                    setAxisLabelsAndRanges(graph, xLabel, yLabels[i], xLimits, 
                                           (suffixes[i] == "ALL") ? std::make_pair(-0.1, 0.6) :
                                           (suffixes[i] == "doubleratio") ? std::make_pair(-0.02, 0.3) :
                                           std::make_pair(-0.1, 0.1));
                }

                graph->Draw((datasetIndex == 0) ? "AP" : "P SAME");
            }
        }

        // Add the dashed gray line at y = 0
        TLine* line = new TLine(xLimits.first, 0, xLimits.second, 0);
        line->SetLineColor(kGray+2);
        line->SetLineStyle(7);  // Dashed line
        line->Draw("same");

        // Add the legend if we're at the last subplot
        if (i == 5) {  // Adjust the position of the legend
            TLegend* legend = new TLegend(0.1, 0.7, 0.48, 0.9);
            legend->SetTextFont(42);
            legend->SetFillColor(0);
            legend->SetBorderSize(1);

            for (size_t j = 0; j < legendEntries.size(); ++j) {
                legend->AddEntry(graphs[j], legendEntries[j].c_str(), "P");
            }

            legend->Draw();
        }
    }

    gSystem->Exec("mkdir -p output/rho0_plots");
    c->SaveAs(outputFileName.c_str());
    // delete c;
}


int main(int argc, char** argv) {

	std::string asymmetryFile = argv[1];
    // Load asymmetry data
    std::map<std::string, std::vector<std::vector<double>>> asymmetryData = readAsymmetries(asymmetryFile);

    // // Output the contents of asymmetryData
    // for (const auto& entry : asymmetryData) {
    //     std::cout << "Key: " << entry.first << std::endl;
    //     const auto& values = entry.second;

    //     for (const auto& row : values) {
    //         std::cout << "    x: " << row[0] 
    //                   << ", y: " << row[1] 
    //                   << ", error: " << row[2] 
    //                   << std::endl;
    //     }
    // }

    // Define the 8 prefixes that correspond to the different datasets
    std::vector<std::string> prefixes = {
        "epiplus", 
        "epipluspiminus", 
        "epipluspiminus_rho0_free", 
        "eppiplus", 
        "eppiplus_rho0_free", 
        "eppipluspiminus", 
        "eppipluspiminus_rho0_free_A", 
        "eppipluspiminus_rho0_free_B"
    };

    // Define the legend entries
    std::vector<std::string> legendEntries = {
        "e#pi^{+}X, M_{x (e' #pi^{+})} > 1.5 GeV",
        "e#pi^{+}#pi^{-}, M_{x (e'#pi^{+}} > 1.5 GeV",
        "e#pi^{+}#pi^{-}, M_{x (e'#pi^{+}} > 1.5 GeV, M_{x (e'#pi^{+}#pi^{-})} > 1.05 GeV",
        "ep#pi^{+}X, M_{x (e'#pi^{+}} > 1.5 GeV",
        "ep#pi^{+}X, M_{x (e'#pi^{+}} > 1.5 GeV, M_{x (e'p} > 0.95 GeV",
        "ep#pi^{+}#pi^{-}X, M_{x (e' #pi^{+})} > 1.5 GeV",
        "ep#pi^{+}#pi^{-}X, M_{x (e' #pi^{+})} > 1.5 GeV, M_{x (e'#pi^{+}#pi^{-})} > 1.05 GeV",
        "ep#pi^{+}#pi^{-}X, M_{x (e' #pi^{+})} > 1.5 GeV, M_{x (e'p} > 0.95 GeV"
    };

    // Call plotDependence
    plotDependence(asymmetryData, prefixes, "P_{T} (GeV)", {0.0, 1.1}, "output/rho0_plots/PT_dependence_plots.png", legendEntries);

    return 0;
}