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
#include <TPad.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TText.h>

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
    const int color
) {
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
    const std::pair<double, double> &yLimits
) {
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

// Function to partition the canvas
void CanvasPartition(TCanvas *C, const Int_t Nx = 2, const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05) {
    if (!C) return;

    Float_t vSpacing = 0.0;
    Float_t vStep  = (1. - bMargin - tMargin - (Ny - 1) * vSpacing) / Ny;
    Float_t hSpacing = 0.0;
    Float_t hStep  = (1. - lMargin - rMargin - (Nx - 1) * hSpacing) / Nx;

    Float_t vposd, vposu, vmard, vmaru, vfactor;
    Float_t hposl, hposr, hmarl, hmarr, hfactor;

    for (Int_t i = 0; i < Nx; i++) {
        if (i == 0) {
            hposl = 0.0;
            hposr = lMargin + hStep;
            hfactor = hposr - hposl;
            hmarl = lMargin / hfactor;
            hmarr = 0.0;
        } else if (i == Nx - 1) {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep + rMargin;
            hfactor = hposr - hposl;
            hmarl = 0.0;
            hmarr = rMargin / hfactor;
        } else {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep;
            hfactor = hposr - hposl;
            hmarl = 0.0;
            hmarr = 0.0;
        }

        for (Int_t j = 0; j < Ny; j++) {
            if (j == 0) {
                vposd = 0.0;
                vposu = bMargin + vStep;
                vfactor = vposu - vposd;
                vmard = bMargin / vfactor;
                vmaru = 0.0;
            } else if (j == Ny - 1) {
                vposd = vposu + vSpacing;
                vposu = vposd + vStep + tMargin;
                vfactor = vposu - vposd;
                vmard = 0.0;
                vmaru = tMargin / vfactor;
            } else {
                vposd = vposu + vSpacing;
                vposu = vposd + vStep;
                vfactor = vposu - vposd;
                vmard = 0.0;
                vmaru = 0.0;
            }

            C->cd(0);
            auto name = TString::Format("pad_%d_%d", i, j);
            auto pad = (TPad*) C->FindObject(name.Data());
            if (pad) delete pad;
            pad = new TPad(name.Data(), "", hposl, vposd, hposr, vposu);
            pad->SetLeftMargin(hmarl);
            pad->SetRightMargin(hmarr);
            pad->SetBottomMargin(vmard);
            pad->SetTopMargin(vmaru);

            pad->SetFrameBorderMode(0);
            pad->SetBorderMode(0);
            pad->SetBorderSize(0);

            pad->Draw();
        }
    }
}

void plotQ2yz_pT(const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData, const std::string &outputFileName) {
    const Int_t Nx = 5;
    const Int_t Ny = 4;

    TCanvas *c = new TCanvas("c", "Q2, y, z Dependence", 1600, 1600);

    // Partition the canvas
    CanvasPartition(c, Nx, Ny, 0.05, 0.05, 0.05, 0.05);

    std::vector<std::string> Q2yPrefixes = {"Q2y1", "Q2y5", "Q2y9", "Q2y13", "Q2y16", 
                                            "Q2y2", "Q2y6", "Q2y10", "Q2y14", "Q2y17", 
                                            "Q2y3", "Q2y7", "Q2y11", "Q2y15", "", 
                                            "Q2y4", "Q2y8", "Q2y12", "", ""};
    std::vector<int> colors = {kBlack, kRed, kGreen + 2, kBlue, kMagenta};

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int index = i + j * Nx;
            if (Q2yPrefixes[index].empty()) continue;

            c->cd(index + 1);
            gPad->SetFillStyle(4000);
            gPad->SetFrameFillStyle(4000);

            bool dataFound = false;

            for (int zIndex = 1; zIndex <= 5; ++zIndex) {
                std::string key = Q2yPrefixes[index] + "z" + std::to_string(zIndex) + "chi2FitsALUsinphi";
                auto it = asymmetryData.find(key);

                if (it != asymmetryData.end()) {
                    dataFound = true;
                    const auto &data = it->second;

                    std::vector<double> x, y, yErr;
                    for (const auto &entry : data) {
                        x.push_back(entry[0]);
                        y.push_back(entry[1]);
                        yErr.push_back(entry[2]);
                    }

                    TGraphErrors *graph = createTGraphErrors(x, y, yErr, 20, 0.8, colors[zIndex - 1]);
                    setAxisLabelsAndRanges(graph, "P_{T} (GeV)", "F_{LU}^{sin#phi}/F_{UU}", {0.0, 1.0}, {-0.1, 0.1});

                    // Hide X and Y axis labels for inner plots
                    if (i != 0) { // Not the first column
                        graph->GetYaxis()->SetLabelOffset(999);
                        graph->GetYaxis()->SetTitleOffset(999);
                    }
                    if (j != Ny - 1) { // Not the last row
                        graph->GetXaxis()->SetLabelOffset(999);
                        graph->GetXaxis()->SetTitleOffset(999);
                    }

                    graph->Draw(zIndex == 1 ? "AP" : "P SAME");
                }
            }

            if (!dataFound) {
                std::cout << "Warning: No data found for prefix " << Q2yPrefixes[index] << std::endl;
            }
        }
    }

    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs(outputFileName.c_str());

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

    // // Call the plotting function for different dependencies
    // plotDependence(asymmetryData, "x", "x_{B}", {0.06, 0.6}, "output/epX_plots/x_dependence_plots.png");
    // plotDependence(asymmetryData, "PT", "P_{T} (GeV)", {0.0, 1.0}, "output/epX_plots/PT_dependence_plots.png");
    // plotDependence(asymmetryData, "xF", "x_{F}", {-0.8, 0.6}, "output/epX_plots/xF_dependence_plots.png");

    // // Plot PT and xF dependence comparison
    // plotComparison(asymmetryData, "output/epX_plots/PT_xF_dependence_comparison.png");

    // Plot Q2-y-z dependence
    plotQ2yz_pT(asymmetryData, "output/epX_plots/Q2yz_dependence_plots.png");

    return 0;
}
