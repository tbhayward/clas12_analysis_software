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


// Function to read arrays directly from the kinematic file
std::map<std::string, std::vector<std::vector<double>>> readKinematics(const std::string &filename) {
    std::map<std::string, std::vector<std::vector<double>>> kinematicData;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::string::size_type pos = line.find("=");
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos - 1);
            std::string dataStr = line.substr(pos + 1);

            // Remove unnecessary characters from the data string
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), '{'), dataStr.end());
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), '}'), dataStr.end());
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), ';'), dataStr.end());

            // Split dataStr into individual numbers
            std::vector<std::vector<double>> values;
            std::stringstream ss(dataStr);
            std::string num;
            std::vector<double> tempVec;
            while (std::getline(ss, num, ',')) {
                if (!num.empty()) {
                    tempVec.push_back(std::stod(num));
                }
                if (tempVec.size() == 9) { // Assuming each kinematic entry has 9 values
                    values.push_back(tempVec);
                    tempVec.clear();
                }
            }

            kinematicData[key] = values;
        }
    }

    return kinematicData;
}

// Function to read asymmetry fits directly from the file
std::map<std::string, std::vector<std::vector<double>>> readAsymmetries(const std::string &filename) {
    std::map<std::string, std::vector<std::vector<double>>> asymmetryData;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::string::size_type pos = line.find("=");
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos - 1);
            std::string dataStr = line.substr(pos + 1);

            // Remove unnecessary characters from the data string
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), '{'), dataStr.end());
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), '}'), dataStr.end());
            dataStr.erase(std::remove(dataStr.begin(), dataStr.end(), ';'), dataStr.end());

            // Split dataStr into individual numbers
            std::vector<std::vector<double>> values;
            std::stringstream ss(dataStr);
            std::string num;
            std::vector<double> tempVec;
            while (std::getline(ss, num, ',')) {
                if (!num.empty()) {
                    tempVec.push_back(std::stod(num));
                }
                if (tempVec.size() == 3) { // Assuming each asymmetry entry has 3 values
                    values.push_back(tempVec);
                    tempVec.clear();
                }
            }

            asymmetryData[key] = values;
        }
    }

    return asymmetryData;
}

// Function to print the data for verification
void printData(const std::map<std::string, std::vector<std::vector<double>>> &data) {
    for (const auto &entry : data) {
        std::cout << entry.first << " = {\n";
        for (const auto &vec : entry.second) {
            std::cout << "  {";
            for (size_t i = 0; i < vec.size(); ++i) {
                std::cout << vec[i];
                if (i < vec.size() - 1) std::cout << ", ";
            }
            std::cout << "},\n";
        }
        std::cout << "};\n";
    }
}

void plotDependence(
    const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData,
    const std::string &prefix, 
    const std::string &xLabel, 
    const std::pair<double, double> &xLimits, 
    const std::string &outputFileName
) {
    // Create a 2x3 canvas
    TCanvas *c = new TCanvas("c", "Dependence Plots", 1200, 800);
    c->Divide(3, 2);

    // Define the suffixes for the asymmetries we want to plot
    std::vector<std::string> suffixes = {"ALUsinphi", "AULsinphi", "AULsin2phi", "ALL", "ALLcosphi"};
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin(2#phi)}/F_{UU}",
        "F_{LL}/F_{UU}",
        "F_{LL}^{cos#phi}/F_{UU}"
    };

    // Plot each asymmetry in its respective subplot
    for (size_t i = 0; i < suffixes.size(); ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.18);  // Increase left margin for Y-axis label
        gPad->SetBottomMargin(0.15);  // Increase bottom margin for X-axis label

        // Construct the key by combining prefix and suffix
        std::string key = prefix + "chi2Fits" + suffixes[i];

        // Check if the key exists in the map
        auto it = asymmetryData.find(key);
        if (it != asymmetryData.end()) {
            const auto &data = it->second;

            // Create vectors to hold x, y, y error (statistical), and combined uncertainties
            std::vector<double> x, y, yStatErr, yCombErr;
            for (size_t j = 0; j < data.size(); ++j) {
                x.push_back(data[j][0]);
                y.push_back(data[j][1]);
                yStatErr.push_back(data[j][2]);

                // Determine systematic uncertainty based on the suffix
                double sysUncertainty = 0.0;
                if (suffixes[i] == "ALUsinphi") {
                    sysUncertainty = y[j] * 0.029;  // 2.9% of the current y-value
                } else if (suffixes[i] == "AULsinphi" || suffixes[i] == "AULsin2phi") {
                    sysUncertainty = y[j] * 0.101;  // 10.1% of the current y-value
                } else if (suffixes[i] == "ALL" || suffixes[i] == "ALLcosphi") {
                    sysUncertainty = std::sqrt(
                        std::pow(y[j] * 0.029, 2) +  // 2.9% of the current y-value
                        std::pow(y[j] * 0.101, 2)   // 10.1% of the current y-value
                    );
                }

                // Combine statistical and systematic uncertainties in quadrature
                yCombErr.push_back(std::sqrt(std::pow(yStatErr[j], 2) + std::pow(sysUncertainty, 2)));
            }

            // Create TGraphErrors for the combined (stat + sys) uncertainties
            TGraphErrors *graphComb = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, yCombErr.data());
            graphComb->SetMarkerStyle(20);  // Circle points
            graphComb->SetMarkerSize(0.8);  // Smaller marker size
            graphComb->SetMarkerColor(kRed-7);  // Light red color
            graphComb->SetLineColor(kRed-7);  // Light red color

            // Create TGraphErrors for the statistical uncertainties
            TGraphErrors *graphStat = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, yStatErr.data());
            graphStat->SetMarkerStyle(20);  // Circle points
            graphStat->SetMarkerSize(0.8);  // Smaller marker size
            graphStat->SetMarkerColor(kBlack);
            graphStat->SetLineColor(kBlack);

            // Set x-axis and y-axis ranges based on the statistical uncertainties graph
            graphStat->GetXaxis()->SetLimits(xLimits.first, xLimits.second);
            graphStat->GetXaxis()->SetRangeUser(xLimits.first, xLimits.second);
            if (suffixes[i] == "ALL") {
                graphStat->GetYaxis()->SetRangeUser(-0.1, 0.6);
            } else {
                graphStat->GetYaxis()->SetRangeUser(-0.15, 0.15);
            }

            // Draw combined uncertainties (statistical + systematic) first
            graphComb->Draw("AP");

            // Redraw the statistical uncertainties on top
            graphStat->Draw("P SAME");

            // Draw a faint dashed gray horizontal line at y=0
            TLine *line = new TLine(graphStat->GetXaxis()->GetXmin(), 0, graphStat->GetXaxis()->GetXmax(), 0);
            line->SetLineColor(kGray+2);
            line->SetLineStyle(7);  // Dashed line
            line->Draw();
        }
    }

    // Save the canvas as a PNG file
    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs(outputFileName.c_str());

    // Clean up
    delete c;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <asymmetries.txt> <kinematicPlots.txt>\n";
        return 1;
    }

    std::string asymmetryFile = argv[1];
    std::string kinematicFile = argv[2];

    // Read the asymmetry data from the file
    std::map<std::string, std::vector<std::vector<double>>> asymmetryData = readAsymmetries(asymmetryFile);

    // Read the kinematic data from the file
    std::map<std::string, std::vector<std::vector<double>>> kinematicData = readKinematics(kinematicFile);

    // // Print out the parsed data
    // std::cout << "Asymmetry Data:\n";
    // printData(asymmetryData);

    // std::cout << "\nKinematic Data:\n";
    // printData(kinematicData);

    // Call the plotting function for different dependencies
    plotDependence(asymmetryData, "x", "x_{B}", {0.06, 0.6}, "output/epX_plots/x_dependence_plots.png");
    plotDependence(asymmetryData, "PT", "P_{T} (GeV)", {0.0, 1.0}, "output/epX_plots/PT_dependence_plots.png");
    plotDependence(asymmetryData, "xF", "x_{F}", {-0.8, 0.8}, "output/epX_plots/xF_dependence_plots.png");

    return 0;
}
