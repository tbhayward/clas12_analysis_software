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
#include <TLegend.h>
#include <TLine.h>
#include <TLegendEntry.h>
#include <TAxis.h>
#include <TSystem.h>
#include <vector>
#include <string>
#include <map>
#include <TLatex.h>
#include <TText.h>  // Add this include for TText

// Global variables for systematic uncertainties
const double LU_SYS_UNCERTAINTY = 0.029;
const double UL_SYS_UNCERTAINTY = 0.101;
const double LL_SYS_UNCERTAINTY_LU_COMPONENT = LU_SYS_UNCERTAINTY;
const double LL_SYS_UNCERTAINTY_UL_COMPONENT = UL_SYS_UNCERTAINTY;


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

TGraphErrors* createTGraphErrors(
    const std::vector<double> &x, 
    const std::vector<double> &y, 
    const std::vector<double> &yErr,
    const int markerStyle, 
    const double markerSize, 
    const int color) {
    TGraphErrors *graph = new TGraphErrors(x.size(), x.data(), y.data(), nullptr, yErr.data());
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetTitle("");
    return graph;
}

void setAxisLabelsAndRanges(
    TGraphErrors *graph, 
    const std::string &xLabel, 
    const std::string &yLabel, 
    const std::pair<double, double> &xLimits, 
    const std::pair<double, double> &yLimits) {
    graph->GetXaxis()->SetTitle(xLabel.c_str());
    graph->GetYaxis()->SetTitle(yLabel.c_str());
    graph->GetXaxis()->SetLimits(xLimits.first, xLimits.second);
    graph->GetXaxis()->SetRangeUser(xLimits.first, xLimits.second);
    graph->GetYaxis()->SetRangeUser(yLimits.first, yLimits.second);

    // Adjust axis label and title sizes
    graph->GetXaxis()->SetLabelSize(0.04);
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetLabelSize(0.04);
    graph->GetYaxis()->SetTitleSize(0.045);
}

void setCustomAxisLabelsAndRanges(
    TGraphErrors *graph, 
    const std::string &xLabel, 
    const std::string &yLabel, 
    const std::pair<double, double> &xLimits, 
    const std::pair<double, double> &yLimits) {
    
    graph->GetXaxis()->SetTitle(xLabel.c_str());
    graph->GetYaxis()->SetTitle(yLabel.c_str());
    graph->GetXaxis()->SetLimits(xLimits.first, xLimits.second);
    graph->GetXaxis()->SetRangeUser(xLimits.first, xLimits.second);
    
    // Set custom y-axis range
    graph->GetYaxis()->SetRangeUser(yLimits.first, yLimits.second);

    // Adjust axis label and title sizes
    graph->GetXaxis()->SetLabelSize(0.04);
    graph->GetXaxis()->SetTitleSize(0.045);
    graph->GetYaxis()->SetLabelSize(0.04);
    graph->GetYaxis()->SetTitleSize(0.045);
}

double computeSystematicUncertainty(const std::string &suffix, double yValue) {
    if (suffix == "ALUsinphi") {
        return yValue * LU_SYS_UNCERTAINTY;  // LU systematic uncertainty
    } else if (suffix == "AULsinphi" || suffix == "AULsin2phi") {
        return yValue * UL_SYS_UNCERTAINTY;  // UL systematic uncertainty
    } else if (suffix == "ALL" || suffix == "ALLcosphi") {
        return std::sqrt(
            std::pow(yValue * LL_SYS_UNCERTAINTY_LU_COMPONENT, 2) +  // LU component
            std::pow(yValue * LL_SYS_UNCERTAINTY_UL_COMPONENT, 2)    // UL component
        );
    }
    return 0.0;
}

void plotDependence(
    const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData,
    const std::string &prefix, 
    const std::string &xLabel, 
    const std::pair<double, double> &xLimits, 
    const std::string &outputFileName) {
    TCanvas *c = new TCanvas("c", "Dependence Plots", 1200, 800);
    c->Divide(3, 2);

    std::vector<std::string> suffixes = {"ALUsinphi", "AULsinphi", "AULsin2phi", "AULoffset", "ALL", "ALLcosphi"};
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin(2#phi)}/F_{UU}",
        "A_{UL} offset",
        "F_{LL}/F_{UU}",
        "F_{LL}^{cos#phi}/F_{UU}"
    };

    for (size_t i = 0; i < suffixes.size(); ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.15);

        std::string key = prefix + "chi2Fits" + suffixes[i];
        auto it = asymmetryData.find(key);
        if (it != asymmetryData.end()) {
            const auto &data = it->second;

            std::vector<double> x, y, yStatErr, yCombErr;
            for (const auto &entry : data) {
                x.push_back(entry[0]);
                y.push_back(entry[1]);
                yStatErr.push_back(entry[2]);

                double sysUncertainty = computeSystematicUncertainty(suffixes[i], y.back());
                yCombErr.push_back(std::sqrt(std::pow(yStatErr.back(), 2) + std::pow(sysUncertainty, 2)));
            }

            TGraphErrors *graphComb = nullptr;
            if (suffixes[i] != "AULoffset") {
                graphComb = createTGraphErrors(x, y, yCombErr, 20, 0.8, kRed-7);
                setAxisLabelsAndRanges(graphComb, xLabel, yLabels[i], xLimits, (suffixes[i] == "ALL") ? std::make_pair(-0.1, 0.6) : std::make_pair(-0.15, 0.15));
                graphComb->Draw("AP");
            }

            TGraphErrors *graphStat = createTGraphErrors(x, y, yStatErr, 20, 0.8, kBlack);
            setAxisLabelsAndRanges(graphStat, xLabel, yLabels[i], xLimits, (suffixes[i] == "AULoffset") ? std::make_pair(-0.2, 0.2) : (suffixes[i] == "ALL") ? std::make_pair(-0.1, 0.6) : std::make_pair(-0.15, 0.15));

            if (suffixes[i] != "AULoffset") {
                graphStat->Draw("P SAME");
            } else {
                graphStat->Draw("AP");
            }

            TLine *line = new TLine(xLimits.first, 0, xLimits.second, 0);
            line->SetLineColor(kGray+2);
            line->SetLineStyle(7);
            line->Draw();
        }
    }

    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs(outputFileName.c_str());
    delete c;
}

void plotComparison(
    const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData,
    const std::string &outputFileName) {
    // Create a 1x2 canvas
    TCanvas *c = new TCanvas("c", "PT and xF Dependence Comparison", 1200, 600);
    c->Divide(2, 1); // 1 row, 2 columns

    // Define the H2 PT dependence data
    std::vector<std::vector<double>> H2DataPT = {
        {0.067870, -0.004157, 0.001025}, {0.155431, -0.007453, 0.000539},
        {0.251797, -0.011001, 0.000442}, {0.349481, -0.016850, 0.000421},
        {0.448136, -0.024485, 0.000449}, {0.547082, -0.033942, 0.000525},
        {0.646201, -0.047847, 0.000661}, {0.745191, -0.058604, 0.000889},
        {0.844268, -0.073812, 0.001278}, {0.943535, -0.074605, 0.001973},
        {1.042894, -0.086452, 0.003307}, {1.142434, -0.074708, 0.006179},
        {1.242379, -0.066435, 0.013280}
    };

    // Define the H2 xF dependence data
    std::vector<std::vector<double>> H2DataXF = {
        {-0.75, -0.002333, 0.007808}, {-0.65, -0.001810, 0.001441},
        {-0.55, -0.017523, 0.000810}, {-0.45, -0.024326, 0.000588},
        {-0.35, -0.027394, 0.000502}, {-0.25, -0.026728, 0.000479},
        {-0.15, -0.024592, 0.000488}, {-0.05, -0.018854, 0.000524},
        {0.05, -0.012483, 0.000597}, {0.15, -0.002537, 0.000723},
        {0.25, 0.009889, 0.000929}, {0.35, 0.022082, 0.001272},
        {0.45, 0.030298, 0.001874}, {0.55, 0.041041, 0.003080},
        {0.65, 0.048530, 0.006303}, {0.75, 0.102060, 0.026785}
    };

    // Define the keys for the NH3 data
    std::string suffix = "ALUsinphi";
    std::string keyPT = "PTchi2Fits" + suffix;
    std::string keyXF = "xFchi2Fits" + suffix;

    // Check if the keys exist in the map for NH3 data
    auto itPT = asymmetryData.find(keyPT);
    auto itXF = asymmetryData.find(keyXF);

    if (itPT != asymmetryData.end() && itXF != asymmetryData.end()) {
        const auto &NH3DataPT = itPT->second;
        const auto &NH3DataXF = itXF->second;

        // Create vectors for NH3 PT data
        std::vector<double> xNH3PT, yNH3PT, yErrNH3PT;
        for (const auto &entry : NH3DataPT) {
            xNH3PT.push_back(entry[0]);
            yNH3PT.push_back(entry[1]);
            yErrNH3PT.push_back(entry[2]);
        }

        // Create vectors for NH3 xF data
        std::vector<double> xNH3XF, yNH3XF, yErrNH3XF;
        for (const auto &entry : NH3DataXF) {
            xNH3XF.push_back(entry[0]);
            yNH3XF.push_back(entry[1]);
            yErrNH3XF.push_back(entry[2]);
        }

        // Create vectors for H2 PT data
        std::vector<double> xH2PT, yH2PT, yErrH2PT;
        for (const auto &entry : H2DataPT) {
            xH2PT.push_back(entry[0]);
            yH2PT.push_back(entry[1]);
            yErrH2PT.push_back(entry[2]);
        }

        // Create vectors for H2 xF data
        std::vector<double> xH2XF, yH2XF, yErrH2XF;
        for (const auto &entry : H2DataXF) {
            xH2XF.push_back(entry[0]);
            yH2XF.push_back(entry[1]);
            yErrH2XF.push_back(entry[2]);
        }

        // Create the NH3 PT TGraphErrors
        TGraphErrors *graphNH3PT = createTGraphErrors(xNH3PT, yNH3PT, yErrNH3PT, 20, 0.8, kRed);
        setAxisLabelsAndRanges(graphNH3PT, "P_{T} (GeV)", "F_{LU}^{sin#phi}/F_{UU}", {0.0, 1.2}, {-0.1, 0.1});

        // Create the H2 PT TGraphErrors
        TGraphErrors *graphH2PT = createTGraphErrors(xH2PT, yH2PT, yErrH2PT, 21, 0.8, kBlue);

        // Create the NH3 xF TGraphErrors
        TGraphErrors *graphNH3XF = createTGraphErrors(xNH3XF, yNH3XF, yErrNH3XF, 20, 0.8, kRed);
        setAxisLabelsAndRanges(graphNH3XF, "x_{F}", "F_{LU}^{sin#phi}/F_{UU}", {-0.8, 0.8}, {-0.1, 0.1});

        // Create the H2 xF TGraphErrors
        TGraphErrors *graphH2XF = createTGraphErrors(xH2XF, yH2XF, yErrH2XF, 21, 0.8, kBlue);

        // Plot PT dependence in the first pad
        c->cd(1);
        graphNH3PT->Draw("AP");
        graphH2PT->Draw("P SAME");

        // Add horizontal line at y=0
        TLine *linePT = new TLine(0, 0, 1.2, 0);
        linePT->SetLineColor(kGray + 2);
        linePT->SetLineStyle(2);
        linePT->Draw("same");

        // Add legend to the first pad
        TLegend *legendPT = new TLegend(0.7, 0.8, 0.9, 0.9);
        legendPT->AddEntry(graphNH3PT, "NH_{3}", "p");
        legendPT->AddEntry(graphH2PT, "H_{2}", "p");
        legendPT->Draw();

        // Plot xF dependence in the second pad
        c->cd(2);
        graphNH3XF->Draw("AP");
        graphH2XF->Draw("P SAME");

        // Add horizontal line at y=0
        TLine *lineXF = new TLine(-0.8, 0, 0.8, 0);
        lineXF->SetLineColor(kGray + 2);
        lineXF->SetLineStyle(2);
        lineXF->Draw("same");

        // Add legend to the second pad
        TLegend *legendXF = new TLegend(0.7, 0.8, 0.9, 0.9);
        legendXF->AddEntry(graphNH3XF, "NH_{3}", "p");
        legendXF->AddEntry(graphH2XF, "H_{2}", "p");
        legendXF->Draw();

        // Save the canvas as a PNG file
        gSystem->Exec("mkdir -p output/epX_plots");
        c->SaveAs(outputFileName.c_str());

        // Clean up
        delete c;
        delete legendPT;
        delete legendXF;
        delete linePT;
        delete lineXF;
    } else {
        if (itPT == asymmetryData.end()) std::cerr << "Error: No NH3 PT data found for key " << keyPT << "\n";
        if (itXF == asymmetryData.end()) std::cerr << "Error: No NH3 xF data found for key " << keyXF << "\n";
    }
}

TCanvas *setupCanvas(int width, int height, int cols, int rows) {
    // Create a canvas with the specified width and height
    TCanvas *c = new TCanvas("c", "Q2-y-z Dependence", width, height);
    
    // Divide the canvas into a grid of cols x rows sub-pads without extra row for titles
    c->Divide(cols, rows, 0, 0);  
    
    return c;
}

// Function to create and draw data plots with titles above the plot area
void drawDataPlotWithTitle(TGraphErrors* graph, int q2Index, int row, bool firstGraphDrawn, const std::string& title = "") {
    if (!firstGraphDrawn) {
        // Adjust the pad's top margin to create space for the title
        if (row == 0) {
            gPad->SetTopMargin(0.1);  // Increase top margin for the top row
        }

        // Draw the plot
        setAxisLabelsAndRanges(graph, "P_{T} (GeV)", "F_{LU}^{sin#phi}/F_{UU}", {0.1, 0.9}, {-0.09, 0.09});
        graph->GetXaxis()->SetLabelFont(42); // 42 is bold Helvetica
        graph->GetYaxis()->SetLabelFont(42);
        graph->GetXaxis()->SetTitleFont(42);
        graph->GetYaxis()->SetTitleFont(42);
        graph->GetXaxis()->SetLabelSize(0.05); 
        graph->GetYaxis()->SetLabelSize(0.05); 
        graph->GetXaxis()->SetTitleSize(0.06); 
        graph->GetYaxis()->SetTitleSize(0.06); 
        graph->GetXaxis()->SetNdivisions(505);

        // Hide Y-axis labels for non-leftmost plots
        if (q2Index != 0) {
            graph->GetYaxis()->SetLabelOffset(999);
            graph->GetYaxis()->SetTitleOffset(999);
        }

        // Hide X-axis labels for non-bottom row plots and specific bottom right subplots
        if (row != 3 || (row == 3 && (q2Index == 3 || q2Index == 4))) {
            graph->GetXaxis()->SetLabelOffset(999);
            graph->GetXaxis()->SetTitleOffset(999);
        }

        graph->Draw("AP");

        // Add title above the plot
        if (!title.empty() && row == 0) {
            TLatex latex;
            latex.SetNDC();  // Use Normalized Device Coordinates
            latex.SetTextAlign(22);  // Center alignment
            latex.SetTextSize(0.06); // Adjust text size as needed
            latex.DrawLatex(0.5, 0.95, title.c_str());  // Positioned above the plot inside the pad
        }
    } else {
        graph->Draw("P SAME");
    }
}

// Function to add the legend
void addLegend(std::vector<TGraph*>& sampleGraphs, TCanvas* c) {
    c->cd(20); // Go to the last pad
    TLegend *legend = new TLegend(0.225, 0.225, 0.9, 0.9); // Adjust position and size of the legend box

    // Define the z ranges explicitly
    std::vector<std::string> zRanges = {
        "0.10 < z < 0.25",
        "0.25 < z < 0.35",
        "0.35 < z < 0.45",
        "0.45 < z < 0.55",
        "0.55 < z < 0.75"
    };

    // Add each z range to the legend
    for (size_t zIndex = 0; zIndex < zRanges.size(); ++zIndex) {
        legend->AddEntry(sampleGraphs[zIndex], zRanges[zIndex].c_str(), "P");
    }

    legend->SetTextSize(0.05); // Adjust text size if needed
    legend->SetFillColor(0); // Make background transparent
    legend->SetLineColor(1); // Add border
    legend->Draw();
}

// Function to handle empty placeholders
void drawEmptyPlot(TGraphErrors* dummyGraph, int q2Index, int row, int totalRows) {
    dummyGraph->Draw("AP");
    dummyGraph->GetXaxis()->SetNdivisions(505);

    // Hide Y-axis labels for non-leftmost plots
    if (q2Index != 0) {
        dummyGraph->GetYaxis()->SetLabelOffset(999);
        dummyGraph->GetYaxis()->SetTitleOffset(999);
    }

    // Hide X-axis labels for non-bottom row plots and specific bottom right subplots
    if (row != totalRows - 1 || (row == 3 && (q2Index == 3 || q2Index == 4))) {
        dummyGraph->GetXaxis()->SetLabelOffset(999);
        dummyGraph->GetXaxis()->SetTitleOffset(999);
    }
}

// Function to filter out data points with large error bars
std::vector<std::vector<double>> filterDataByError(const std::vector<std::vector<double>>& data, double maxError) {
    std::vector<std::vector<double>> filteredData;

    for (const auto& entry : data) {
        double yError = entry[2]; // Assuming the error is stored in the third element of each entry
        if (yError <= maxError) {
            filteredData.push_back(entry); // Include only if error is below threshold
        }
    }

    return filteredData;
}

void addCanvasSideLabels(TCanvas* c, const std::vector<std::string>& y_ranges) {
    c->cd(); // Switch to the entire canvas
    TLatex latex;
    latex.SetNDC();  // Use Normalized Device Coordinates
    latex.SetTextAlign(22);  // Centered alignment
    latex.SetTextSize(0.015);  // Adjust text size
    latex.SetTextAngle(270);  // Rotate text vertically

    // Position and draw each y-range label
    for (size_t i = 0; i < y_ranges.size(); ++i) {
        double yPos = 1.0 - (i + 0.5) * (1.0 / y_ranges.size());
        latex.DrawLatex(0.99, yPos, y_ranges[i].c_str());  // Position text to the right of the plots
    }
}

void plotQ2yz_pT(const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData) {
    std::vector<std::string> fitTypes = {"ALUsinphi", "AULoffset", "AULsinphi", "AULsin2phi", "ALL", "ALLcosphi"};
    std::vector<std::string> yLabels = {"F_{LU}^{sin#phi}/F_{UU}", "A_{UL}^{offset}", "A_{UL}^{sin#phi}", "A_{UL}^{sin2#phi}", "A_{LL}", "A_{LL}^{cos#phi}"};
    std::vector<std::string> outputFiles = {
        "output/epX_plots/Q2yz_pT_ALUsinphi.png",
        "output/epX_plots/Q2yz_pT_AULoffset.png",
        "output/epX_plots/Q2yz_pT_AULsinphi.png",
        "output/epX_plots/Q2yz_pT_AULsin2phi.png",
        "output/epX_plots/Q2yz_pT_ALL.png",
        "output/epX_plots/Q2yz_pT_ALLcosphi.png"
    };

    std::vector<double> maxErrors = {0.0275, 0.01, 0.0275, 0.0275, 0.0275, 0.0275}; 

    std::vector<std::pair<double, double>> yRangesPerPlot = {
        {-0.099, 0.099}, 
        {-0.149, 0.049}, 
        {-0.099, 0.099}, 
        {-0.099, 0.099}, 
        {-0.199, 0.599}, 
        {-0.199, 0.199}  
    };

    TLegend *legend = new TLegend(0.225, 0.225, 0.9, 0.9); 
    std::vector<std::string> zRanges = {
        "0.10 < z < 0.25",
        "0.25 < z < 0.35",
        "0.35 < z < 0.45",
        "0.45 < z < 0.55",
        "0.55 < z < 0.75"
    };
    std::vector<int> colors = {kBlack, kRed, kGreen, kBlue, kMagenta}; 

    for (size_t zIndex = 0; zIndex < zRanges.size(); ++zIndex) {
        TGraph *dummyGraph = new TGraph();
        dummyGraph->SetMarkerColor(colors[zIndex]);
        dummyGraph->SetMarkerStyle(20);
        dummyGraph->SetMarkerSize(1.5);
        legend->AddEntry(dummyGraph, zRanges[zIndex].c_str(), "P");

        TLegendEntry *entry = (TLegendEntry*)legend->GetListOfPrimitives()->Last();
        entry->SetTextColor(colors[zIndex]);
    }
    legend->SetTextSize(0.05); 
    legend->SetFillColor(0);   
    legend->SetLineColor(1);   

    for (size_t fitIndex = 0; fitIndex < fitTypes.size(); ++fitIndex) {
        TCanvas *c = setupCanvas(2400, 1600, 5, 4);

        std::vector<std::vector<std::string>> Q2_prefixes = {
            {"Q2y1", "Q2y5", "Q2y9", "Q2y13", "Q2y16"},
            {"Q2y2", "Q2y6", "Q2y10", "Q2y14", "Q2y17"},
            {"Q2y3", "Q2y7", "Q2y11", "Q2y15", "EMPTY"},
            {"Q2y4", "Q2y8", "Q2y12", "EMPTY", "EMPTY"}
        };
        std::vector<std::string> z_prefixes = {"z1", "z2", "z3", "z4", "z5"};

        std::vector<std::string> topRowTitles = {
            "1.0 < Q^{2} (GeV^{2}) < 2.0",
            "2.0 < Q^{2} (GeV^{2}) < 3.0",
            "3.0 < Q^{2} (GeV^{2}) < 4.0",
            "4.0 < Q^{2} (GeV^{2}) < 5.0",
            "5.0 < Q^{2} (GeV^{2}) < 7.0"
        };

        std::vector<std::string> yRanges = {
            "0.65 < y < 0.75",
            "0.55 < y < 0.65",
            "0.45 < y < 0.55",
            "0.30 < y < 0.45"
        };

        double maxError = maxErrors[fitIndex];
        std::pair<double, double> yRange = yRangesPerPlot[fitIndex];

        for (size_t row = 0; row < Q2_prefixes.size(); ++row) {
            for (size_t q2Index = 0; q2Index < Q2_prefixes[row].size(); ++q2Index) {
                int padIndex = row * 5 + q2Index + 1;
                c->cd(padIndex); 

                if (q2Index != 0) {
                    gPad->SetLeftMargin(0.001);
                } else {
                    gPad->SetLeftMargin(0.18);
                }
                if (row != Q2_prefixes.size() - 1) {
                    gPad->SetBottomMargin(0.001);
                } else {
                    gPad->SetBottomMargin(0.15);
                }

                bool firstGraphDrawn = false;
                bool anyGraphDrawn = false;

                for (size_t zIndex = 0; zIndex < z_prefixes.size(); ++zIndex) {
                    std::string key = Q2_prefixes[row][q2Index] + z_prefixes[zIndex] + "chi2Fits" + fitTypes[fitIndex];
                    auto it = asymmetryData.find(key);

                    if (it == asymmetryData.end()) {
                        continue;
                    }

                    const auto &data = it->second;

                    auto filteredData = filterDataByError(data, maxError);

                    if (filteredData.empty()) continue;

                    std::vector<double> x, y, yErr;

                    for (const auto &entry : filteredData) {
                        x.push_back(entry[0]);
                        y.push_back(entry[1]);
                        yErr.push_back(entry[2]);
                    }

                    TGraphErrors *graph = createTGraphErrors(x, y, yErr, 20, 0.8, colors[zIndex]);

                    // Set the y-axis range directly
                    graph->GetYaxis()->SetRangeUser(yRange.first, yRange.second);

                    std::string title = (row == 0) ? topRowTitles[q2Index] : "";
                    drawDataPlotWithTitle(graph, q2Index, row, firstGraphDrawn, title);
                    firstGraphDrawn = true;
                    anyGraphDrawn = true;
                }

                if (!anyGraphDrawn) {
                    std::vector<double> dummyX = {-9999};
                    std::vector<double> dummyY = {0};
                    std::vector<double> dummyYErr = {0};
                    TGraphErrors *dummyGraph = createTGraphErrors(dummyX, dummyY, dummyYErr, 20, 0.8, kWhite);
                    // setCustomAxisLabelsAndRanges(dummyGraph, "P_{T} (GeV)", yLabels[fitIndex], {0.1, 0.9}, yRange);
                    drawEmptyPlot(dummyGraph, q2Index, row, Q2_prefixes.size());
                }

                if (!(row == 2 && q2Index == 4) && (row != 3 || (q2Index != 3 && q2Index != 4))) {
                    TLine *line = new TLine(0.15, 0.0, 0.95, 0.0);
                    line->SetLineColor(kGray + 2);
                    line->SetLineStyle(7);
                    line->Draw("same");
                }
            }
        }

        addCanvasSideLabels(c, yRanges);

        c->cd(20);
        legend->Draw();

        gSystem->Exec("mkdir -p output/epX_plots");
        c->SaveAs(outputFiles[fitIndex].c_str());
        delete c;
    }

    delete legend;
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

    // // Call the plotting function for different dependencies
    // plotDependence(asymmetryData, "x", "x_{B}", {0.06, 0.6}, "output/epX_plots/x_dependence_plots.png");
    // plotDependence(asymmetryData, "PT", "P_{T} (GeV)", {0.0, 1.0}, "output/epX_plots/PT_dependence_plots.png");
    // plotDependence(asymmetryData, "xF", "x_{F}", {-0.8, 0.6}, "output/epX_plots/xF_dependence_plots.png");

    // // Plot PT and xF dependence comparison
    // plotComparison(asymmetryData, "output/epX_plots/PT_xF_dependence_comparison.png");

    // Plot Q2-y-z dependence
    plotQ2yz_pT(asymmetryData);

    return 0;
}
