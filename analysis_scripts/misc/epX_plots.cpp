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
    const std::string &xLabel, 
    const std::string &outputFileName
) {
    // Create a 2x3 canvas
    TCanvas *c = new TCanvas("c", "Dependence Plots", 1200, 800);
    c->Divide(3, 2);

    // Define the keys for the asymmetries we want to plot
    std::vector<std::string> keys = {"xchi2FitsALUsinphi", "xchi2FitsAULsinphi", "xchi2FitsAULsin2phi", "xchi2FitsALL", "xchi2FitsALLcosphi"};
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin(2#phi)}/F_{UU}",
        "F_{LL}/F_{UU}",
        "F_{LL}^{cos#phi}/F_{UU}"
    };

    // Plot each asymmetry in its respective subplot
    for (size_t i = 0; i < keys.size(); ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.18);  // Increase left margin for Y-axis label
        gPad->SetBottomMargin(0.15);  // Increase bottom margin for X-axis label

        // Check if the key exists in the map
        auto it = asymmetryData.find(keys[i]);
        if (it != asymmetryData.end()) {
            const auto &data = it->second;

            // Create vectors to hold x, y, and y error
            std::vector<double> x, y, yErr;
            for (const auto &entry : data) {
                x.push_back(entry[0]);
                y.push_back(entry[1]);
                yErr.push_back(entry[2]);
            }

            // Create TGraphErrors and plot
            TGraphErrors *graph = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, yErr.data());
            graph->SetTitle("");
            graph->GetXaxis()->SetTitle(xLabel.c_str());
            graph->GetYaxis()->SetTitle(yLabels[i].c_str());

            // Set x-axis and y-axis ranges
            graph->GetXaxis()->SetLimits(0.0, 0.6);
            graph->GetHistogram()->GetXaxis()->SetRangeUser(0.0, 0.6);
            if (keys[i] == "xchi2FitsALL") {
                graph->GetYaxis()->SetRangeUser(-0.1, 0.6);
            } else {
                graph->GetYaxis()->SetRangeUser(-0.15, 0.15);
            }

            // Customize the graph
            graph->SetMarkerStyle(20);  // Circle points
            graph->SetMarkerSize(0.8);  // Smaller marker size
            graph->SetMarkerColor(kBlack);
            graph->SetLineColor(kBlack);

            // Increase the font size of x and y labels
            graph->GetXaxis()->SetLabelSize(0.04);
            graph->GetXaxis()->SetTitleSize(0.05);
            graph->GetYaxis()->SetLabelSize(0.04);
            graph->GetYaxis()->SetTitleSize(0.05);

            graph->Draw("AP");

            // Draw a faint dashed gray horizontal line at y=0
            TLine *line = new TLine(graph->GetXaxis()->GetXmin(), 0, graph->GetXaxis()->GetXmax(), 0);
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
    plotDependence(asymmetryData, "x_{B}", "output/epX_plots/x_dependence_plots.png");
    plotDependence(asymmetryData, "P_{T} (GeV)", "output/epX_plots/PT_dependence_plots.png");
    plotDependence(asymmetryData, "x_{F}", "output/epX_plots/xF_dependence_plots.png");

    return 0;
}
