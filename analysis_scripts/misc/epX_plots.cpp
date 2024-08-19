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
#include <TText.h>  
#include <TF1.h>        
#include <TMarker.h>  
#include <TPaveText.h>

std::vector<std::tuple<int, double, double>> targetPolarizationData = {{{16137, 0.7677, 0.0184}, {16138, 0.7226, 0.0179}, {16144, 0.7484, 0.0217}, {16145, 0.789, 0.0223}, {16146, 0.7063, 0.025}, {16148, 0.7459, 0.0219}, {16156, -0.6508, 0.0519}, {16157, -0.7146, 0.0293}, {16158, -0.7163, 0.0321}, {16164, -0.645, 0.1282}, {16166, -0.779, 0.0375}, {16167, -0.7985, 0.0265}, {16168, -0.779, 0.0788}, {16169, -0.7351, 0.0531}, {16170, -0.674, 0.0211}, {16178, -0.7082, 0.0247}, {16211, 0.8711, 0.0239}, {16214, 0.8096, 0.0231}, {16221, 0.8168, 0.023}, {16222, 0.808, 0.0217}, {16223, 0.7966, 0.0209}, {16224, 0.7067, 0.0682}, {16225, 0.7696, 0.0327}, {16226, 0.893, 0.0365}, {16228, 0.8948, 0.0284}, {16231, -0.7252, 0.0271}, {16232, -0.731, 0.0209}, {16233, -0.7518, 0.0205}, {16234, -0.9052, 0.0398}, {16235, -0.6254, 0.0528}, {16236, -0.6591, 0.0337}, {16238, -0.5733, 0.0488}, {16243, -0.7715, 0.0244}, {16244, -0.7191, 0.0267}, {16245, -0.7029, 0.0196}, {16246, -0.6532, 0.0195}, {16248, -0.6896, 0.0205}, {16249, -0.6618, 0.0394}, {16250, -0.6346, 0.0446}, {16251, -0.8599, 0.0379}, {16252, -0.7241, 0.0331}, {16253, -0.8366, 0.0242}, {16256, -0.6868, 0.0203}, {16257, -0.696, 0.0202}, {16259, -0.6676, 0.0206}, {16260, -0.7266, 0.0254}, {16318, 0.7732, 0.023}, {16320, 0.7917, 0.0216}, {16321, 0.749, 0.0223}, {16322, 0.7914, 0.0285}, {16323, 0.8742, 0.0365}, {16325, 0.677, 0.0643}, {16326, 0.8518, 0.0455}, {16327, 0.7927, 0.0209}, {16328, 0.766, 0.0213}, {16329, 0.7746, 0.0204}, {16330, 0.78, 0.02}, {16331, 0.7265, 0.0293}, {16332, 0.7157, 0.0687}, {16333, 0.7987, 0.0435}, {16335, -0.7002, 0.0214}, {16336, -0.7016, 0.0266}, {16337, -0.6994, 0.0295}, {16338, -0.7127, 0.0222}, {16339, -0.8082, 0.0361}, {16341, -0.6767, 0.0208}, {16343, -0.7024, 0.0213}, {16345, -0.5727, 0.0329}, {16346, -0.7184, 0.0203}, {16348, -0.7732, 0.0211}, {16350, -0.7161, 0.0206}, {16352, -0.7331, 0.0212}, {16353, -0.6895, 0.0199}, {16354, -0.7316, 0.0209}, {16355, -0.7127, 0.0207}, {16356, -0.6927, 0.0208}, {16357, -0.6967, 0.024}, {16709, 0.6617, 0.0334}, {16710, 0.8414, 0.0314}, {16711, 0.7551, 0.0205}, {16712, 0.7663, 0.0204}, {16713, 0.7704, 0.0204}, {16715, 0.7103, 0.0247}, {16716, 0.7885, 0.0268}, {16717, 0.7316, 0.0216}, {16718, 0.6173, 0.0345}, {16719, 0.6784, 0.0205}, {16720, 0.7184, 0.0251}, {16721, -0.4516, 0.068}, {16722, -0.7295, 0.048}, {16723, -0.7011, 0.0219}, {16726, -0.6371, 0.0188}, {16727, -0.6975, 0.0211}, {16728, -0.6279, 0.0195}, {16729, -0.6943, 0.0388}, {16730, -0.6534, 0.0236}, {16731, -0.7173, 0.0276}, {16732, -0.6901, 0.0244}, {16733, -0.7308, 0.0484}, {16734, -0.5724, 0.0968}, {16736, -0.6599, 0.02}, {16738, -0.6721, 0.0239}, {16743, -0.7407, 0.0215}, {16744, -0.6661, 0.0191}, {16746, -0.6751, 0.0196}, {16747, -0.6623, 0.0193}, {16748, -0.6544, 0.0198}, {16749, -0.6525, 0.0438}, {16750, -0.7157, 0.0201}, {16751, -0.6532, 0.0378}, {16752, -0.6773, 0.0301}, {16753, -0.8703, 0.0337}, {16754, -0.6698, 0.0194}, {16755, -0.6767, 0.0195}, {16756, -0.6443, 0.0377}, {16757, -0.6146, 0.0191}, {16759, -0.6446, 0.0187}, {16761, -0.6492, 0.0206}, {16762, -0.6736, 0.0198}, {16763, -0.6451, 0.0244}, {16765, -0.6837, 0.0212}, {16766, -0.6931, 0.0298}, {16767, 0.8475, 0.0223}, {16768, 0.7774, 0.0212}, {16769, 0.7328, 0.0272}, {16770, 0.7803, 0.0219}, {16771, 0.7577, 0.0204}, {16772, 0.7756, 0.022}}};

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

            // Automatically add 'doubleratio' kinematic data if the key is 'ALL'
            if (key.find("ALL") != std::string::npos) {
                std::string doubleRatioKey = key;
                doubleRatioKey.replace(doubleRatioKey.find("ALL"), 3, "doubleratio"); // Replace 'ALL' with 'doubleratio'
                kinematicData[doubleRatioKey] = values; // Copy the same values
            }
        }
    }

    return kinematicData;
}

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

    // Now, calculate the doubleratio fits
    std::string alusKeyPrefix = "ALUsinphi";
    std::string allKeyPrefix = "ALL";
    for (const auto &entry : asymmetryData) {
        // Check if this is an ALUsinphi entry
        if (entry.first.find(alusKeyPrefix) != std::string::npos) {
            std::string baseKey = entry.first.substr(0, entry.first.find(alusKeyPrefix));
            std::string alusKey = baseKey + alusKeyPrefix;
            std::string allKey = baseKey + allKeyPrefix;

            // Ensure both ALUsinphi and ALL exist for this bin
            if (asymmetryData.find(alusKey) != asymmetryData.end() && asymmetryData.find(allKey) != asymmetryData.end()) {
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
                    ratioValue=-ratioValue;
                    doubleratioData.push_back({xValue, ratioValue, error});
                }

                // Store the doubleratio data
                asymmetryData[baseKey + "doubleratio"] = doubleratioData;
            }
        }
    }

    return asymmetryData;
}

// Function to extract and return Q² dependence vectors
std::map<std::string, std::map<std::string, std::vector<std::vector<double>>>>
extractQ2Dependence(
    const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData,
    const std::map<std::string, std::vector<std::vector<double>>> &kinematicData) {

    // Define the fit types we're interested in
    std::vector<std::string> fitTypes = {
        "ALUsinphi", "AULoffset", "AULsinphi", "AULsin2phi", "ALL", "ALLcosphi", "doubleratio"
    };

    // Define the binning structure (prefixes) for each of the cases
    std::map<std::string, std::vector<std::string>> binPrefixes = {
        {"z1pT2y1", {"Q2y1z1", "Q2y5z1", "Q2y9z1", "Q2y13z1", "Q2y16z1"}},
        {"z2pT2y1", {"Q2y1z2", "Q2y5z2", "Q2y9z2", "Q2y13z2", "Q2y16z2"}},
        {"z1pT2y2", {"Q2y2z1", "Q2y6z1", "Q2y10z1", "Q2y14z1", "Q2y17z1"}},
        {"z2pT2y2", {"Q2y2z2", "Q2y6z2", "Q2y10z2", "Q2y14z2", "Q2y17z2"}}
    };

    // Container to store all vectors
    std::map<std::string, std::map<std::string, std::vector<std::vector<double>>>> allVectors;

    // Helper lambda to extract a point
    auto extractPoint = [&](const std::string &prefix, int pointIndex, const std::string &fitType) {
        std::vector<double> result;

        // Get the kinematic data
        if (kinematicData.find(prefix + "Kinematics") != kinematicData.end()) {
            const auto &kinDataVec = kinematicData.at(prefix + "Kinematics");
            double Q2 = kinDataVec[pointIndex][0]; // Q² value
            result.push_back(Q2);
        } else {
            std::cerr << "Kinematic data for " << prefix << " not found.\n";
            return result;
        }

        // Get the asymmetry data
        std::string chi2Fits = prefix + "chi2Fits" + fitType;
        if (asymmetryData.find(chi2Fits) != asymmetryData.end()) {
            const auto &asymDataVec = asymmetryData.at(chi2Fits);
            double asyValue = asymDataVec[pointIndex][1]; // asymmetry value
            double asyErr = asymDataVec[pointIndex][2];   // asymmetry error
            result.push_back(asyValue);
            result.push_back(asyErr);
        } else {
            std::cerr << "Asymmetry data for " << chi2Fits << " not found.\n";
        }

        return result;
    };

    // Extract points for each combination of prefix and fit type
    for (const auto &binPrefix : binPrefixes) {
        const std::string &vectorName = binPrefix.first;
        const std::vector<std::string> &prefixes = binPrefix.second;

        for (const auto &fitType : fitTypes) {
            std::vector<std::vector<double>> dataVector;
            for (const std::string &prefix : prefixes) {
                std::vector<double> point = extractPoint(prefix, 1, fitType); // Use the 2nd point (index 1)
                if (!point.empty()) {
                    dataVector.push_back(point);
                }
            }
            allVectors[vectorName][fitType] = dataVector;
        }
    }

    return allVectors;
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
    const std::string &outputFileName,
    const std::string &extraPrefix = "") {  // Optional extra prefix
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

                double sysUncertainty = 0;  // Optionally calculate systematic uncertainty
                yCombErr.push_back(std::sqrt(std::pow(yStatErr.back(), 2) + std::pow(sysUncertainty, 2)));
            }

            TGraphErrors *graphStat = createTGraphErrors(x, y, yStatErr, 20, 0.8, kBlack);
            setAxisLabelsAndRanges(graphStat, xLabel, yLabels[i], xLimits, 
                                   (suffixes[i] == "AULoffset") ? std::make_pair(-0.2, 0.2) : 
                                   (suffixes[i] == "ALL") ? std::make_pair(-0.1, 0.8) : std::make_pair(-0.08, 0.08));
            graphStat->Draw("AP");

            // Draw the second dataset if the extraPrefix is provided
            if (!extraPrefix.empty()) {
                std::string extraKey = extraPrefix + "chi2Fits" + suffixes[i];
                auto extraIt = asymmetryData.find(extraKey);
                if (extraIt != asymmetryData.end()) {
                    const auto &extraData = extraIt->second;

                    std::vector<double> extraX, extraY, extraYStatErr, extraYCombErr;
                    for (const auto &entry : extraData) {
                        extraX.push_back(entry[0]);
                        extraY.push_back(entry[1]);
                        extraYStatErr.push_back(entry[2]);

                        double extraSysUncertainty = 0;  // Optionally calculate systematic uncertainty
                        extraYCombErr.push_back(std::sqrt(std::pow(extraYStatErr.back(), 2) + std::pow(extraSysUncertainty, 2)));
                    }

                    TGraphErrors *extraGraphStat = createTGraphErrors(extraX, extraY, extraYStatErr, 20, 0.8, kRed);
                    extraGraphStat->Draw("P SAME");
                }

                // Add the text box for labels only if the extraPrefix is provided
                TPaveText *text = new TPaveText(0.65, 0.7, 0.9, 0.9, "NDC");
                if (suffixes[i] == "ALL") {
                    text->SetX1NDC(0.15);  // Move to top left for "F_{LL}/F_{UU}" plot
                    text->SetX2NDC(0.35);
                    text->SetY1NDC(0.75);  // Ensure it's at the top left
                    text->SetY2NDC(0.9);
                }
                text->SetTextAlign(13);
                text->SetBorderSize(1);  // Set border size to 1 for a black border
                text->SetFillColor(0);
                text->AddText("#font[42]{M_{x} > 1.35 GeV}");  // Black text line
                text->AddText("#font[42]{#color[2]{M_{x} > 0 GeV}}");  // Red text line
                text->Draw();
            }
        }
    }

    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs(outputFileName.c_str());
    delete c;
}

void plotRunnumDependence(
    const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData,
    const std::string &prefix, 
    const std::string &xLabel, 
    const std::string &outputFileName) {

    // Create a canvas with 1 row and 2 columns
    TCanvas *c = new TCanvas("c", "Run Number Dependence Plots", 1600, 800);
    c->Divide(2, 1);

    // Define the suffix for FLUsinphi and y-axis label
    std::string suffix = "ALUsinphi";
    std::string yLabel = "F_{LU}^{sin#phi}/F_{UU}";
    std::pair<double, double> yLimits = {-0.06, 0.06};

    // Prepare data vectors for the left and right plots
    std::vector<double> runNumbers, xValues, asymmetries, errors;
    std::vector<double> outlierX, outlierY, outlierErrors; // For outliers

    std::string key = prefix + "chi2Fits" + suffix;
    auto it = asymmetryData.find(key);
    if (it != asymmetryData.end()) {
        const auto &data = it->second;
        for (size_t i = 0; i < data.size(); ++i) {
            runNumbers.push_back(data[i][0]);  // Original x-values (run numbers)
            xValues.push_back(i + 1);          // Sequential run index
            asymmetries.push_back(data[i][1]); // Asymmetry values
            errors.push_back(data[i][2]);      // Asymmetry errors
        }
    }

    // Create a fit function (constant fit)
    TF1 *fitFunc = new TF1("fitFunc", "[0]", xValues.front(), xValues.back());
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineStyle(2);  // Dashed line

    // Fit the data to a constant and retrieve the fit parameters
    TGraphErrors *graphIndex = new TGraphErrors(xValues.size(), xValues.data(), asymmetries.data(), nullptr, errors.data());
    graphIndex->Fit(fitFunc, "Q");  // Silent mode fit
    double mu = fitFunc->GetParameter(0);
    double sigma = fitFunc->GetParError(0);
    double chi2Ndf = fitFunc->GetChisquare() / fitFunc->GetNDF();

    // Draw the left plot: Original run numbers
    c->cd(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);

    TGraphErrors *graphRunnum = createTGraphErrors(runNumbers, asymmetries, errors, 20, 0.8, kBlack);
    setAxisLabelsAndRanges(graphRunnum, xLabel, yLabel, {16135, 16774}, yLimits);
    graphRunnum->Draw("AP");

    // Draw the fitted constant line on the left plot
    fitFunc->Draw("same");

    // Draw the right plot: Run index (1, 2, 3, ...)
    c->cd(2);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);

    // Separate out outliers for the right plot
    for (size_t i = 0; i < asymmetries.size(); ++i) {
        if (std::abs(asymmetries[i] - mu) > 2.5 * errors[i]) {
            // Print outliers
            std::cout << "Outlier found: Run Number " << runNumbers[i] << " " << (std::abs(asymmetries[i] - mu))/errors[i] << std::endl;

            // Store outlier values
            outlierX.push_back(xValues[i]);
            outlierY.push_back(asymmetries[i]);
            outlierErrors.push_back(errors[i]);
        }
    }

    // Create two TGraphErrors: one for regular points, one for outliers
    TGraphErrors *graphRegular = createTGraphErrors(xValues, asymmetries, errors, 20, 0.8, kBlack);
    TGraphErrors *graphOutliers = nullptr;
    if (!outlierX.empty()) {
        graphOutliers = createTGraphErrors(outlierX, outlierY, outlierErrors, 20, 0.8, kRed);
    }

    // Set axis labels and ranges
    setAxisLabelsAndRanges(graphRegular, "run index", yLabel, {0, static_cast<double>(xValues.size()) + 1}, yLimits);
    graphRegular->Draw("AP");

    // Draw the fitted constant line on the right plot
    fitFunc->Draw("same");

    // Draw outliers on the right plot
    if (graphOutliers) {
        graphOutliers->Draw("P SAME");
    }

    // Create and draw a text box in the top right corner with mu, sigma, and chi2/ndf
    for (int i = 1; i <= 2; ++i) {
        c->cd(i);

        TLatex *text = new TLatex();
        text->SetNDC();
        text->SetTextSize(0.0275);
        text->DrawLatex(0.7, 0.85, Form("#mu = %.4g", mu));
        text->DrawLatex(0.7, 0.80, Form("#sigma = %.4g", sigma));
        text->DrawLatex(0.7, 0.75, Form("#chi^{2}/ndf = %.4g", chi2Ndf));
    }

    // Save the canvas to a file
    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs(outputFileName.c_str());

    // Clean up
    delete c;
    delete fitFunc;
}

void plotTargetPolarizationDependence(
    const std::vector<std::tuple<int, double, double>> &targetPolarizationData, 
    const std::string &xLabel, 
    const std::string &outputFileName) {

    // Create a canvas with 1 row and 2 columns
    TCanvas *c = new TCanvas("c", "Target Polarization Dependence Plots", 1600, 800);
    c->Divide(2, 1);

    // Define y-axis label and limits
    std::string yLabel = "Target Polarization";
    std::pair<double, double> yLimits = {-1.0, 1.0};

    // Prepare data vectors for the left and right plots
    std::vector<double> runNumbers, xValues, polarizations, errors;
    std::vector<double> posPolarizations, posErrors, posXValues, posRunNumbers;
    std::vector<double> negPolarizations, negErrors, negXValues, negRunNumbers;

    // Separate the target polarization data into positive and negative values
    for (const auto &entry : targetPolarizationData) {
        int runNumber;
        double polarization, error;
        std::tie(runNumber, polarization, error) = entry;

        runNumbers.push_back(runNumber);
        xValues.push_back(runNumbers.size());  // Sequential run index
        polarizations.push_back(polarization);
        errors.push_back(error);

        if (polarization > 0) {
            posPolarizations.push_back(polarization);
            posErrors.push_back(error);
            posXValues.push_back(xValues.back());
            posRunNumbers.push_back(runNumber);
        } else {
            negPolarizations.push_back(polarization);
            negErrors.push_back(error);
            negXValues.push_back(xValues.back());
            negRunNumbers.push_back(runNumber);
        }
    }

    TF1 *fitFuncPos = new TF1("fitFuncPos", "[0]", 0, 200);
    fitFuncPos->SetLineColor(kRed);
    fitFuncPos->SetLineStyle(2);  // Dashed line

    TF1 *fitFuncNeg = new TF1("fitFuncNeg", "[0]", 0, 200);
    fitFuncNeg->SetLineColor(kBlue);
    fitFuncNeg->SetLineStyle(2);  // Dashed line

    // Perform the fits
    TGraphErrors *graphPos = new TGraphErrors(posXValues.size(), posXValues.data(), posPolarizations.data(), nullptr, posErrors.data());
    TGraphErrors *graphNeg = new TGraphErrors(negXValues.size(), negXValues.data(), negPolarizations.data(), nullptr, negErrors.data());

    graphPos->Fit(fitFuncPos, "Q");  // Silent mode fit
    graphNeg->Fit(fitFuncNeg, "Q");  // Silent mode fit

    double muPos = fitFuncPos->GetParameter(0);
    double sigmaPos = fitFuncPos->GetParError(0);
    double chi2NdfPos = fitFuncPos->GetChisquare() / fitFuncPos->GetNDF();

    double muNeg = fitFuncNeg->GetParameter(0);
    double sigmaNeg = fitFuncNeg->GetParError(0);
    double chi2NdfNeg = fitFuncNeg->GetChisquare() / fitFuncNeg->GetNDF();

    // Draw the left plot: Original run numbers
    c->cd(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);

    TGraphErrors *graphRunnum = createTGraphErrors(runNumbers, polarizations, errors, 20, 0.8, kBlack);
    setAxisLabelsAndRanges(graphRunnum, xLabel, yLabel, {16135, 16774}, yLimits);
    graphRunnum->Draw("AP");

    // Draw the fitted constant lines
    fitFuncPos->Draw("same");
    fitFuncNeg->Draw("same");

    // Draw the right plot: Run index (1, 2, 3, ...)
    c->cd(2);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);

    // Prepare outlier vectors for the right plot
    std::vector<double> outlierPosX, outlierPosY, outlierPosErrors;
    std::vector<double> outlierNegX, outlierNegY, outlierNegErrors;

    for (size_t i = 0; i < polarizations.size(); ++i) {
        double polarization = polarizations[i];
        double error = errors[i];

        if (polarization > 0) {
            if (std::abs(polarization - muPos) > 5.0 * error) {
                outlierPosX.push_back(xValues[i]);
                outlierPosY.push_back(polarization);
                outlierPosErrors.push_back(error);
                std::cout << "Outlier (Positive) found: Run Number " << runNumbers[i] << ", Deviation: " << std::abs(polarization - muPos) / error << std::endl;
            }
        } else {
            if (std::abs(polarization - muNeg) > 5.0 * error) {
                outlierNegX.push_back(xValues[i]);
                outlierNegY.push_back(polarization);
                outlierNegErrors.push_back(error);
                std::cout << "Outlier (Negative) found: Run Number " << runNumbers[i] << ", Deviation: " << std::abs(polarization - muNeg) / error << std::endl;
            }
        }
    }

    // Create TGraphErrors for outliers
    TGraphErrors *graphPosOutliers = nullptr;
    TGraphErrors *graphNegOutliers = nullptr;
    if (!outlierPosX.empty()) {
        graphPosOutliers = createTGraphErrors(outlierPosX, outlierPosY, outlierPosErrors, 20, 0.8, kRed);
    }
    if (!outlierNegX.empty()) {
        graphNegOutliers = createTGraphErrors(outlierNegX, outlierNegY, outlierNegErrors, 20, 0.8, kBlue);
    }

    // Draw the positive and negative polarization data
    TGraphErrors *graphRegularPos = createTGraphErrors(posXValues, posPolarizations, posErrors, 20, 0.8, kBlack);
    TGraphErrors *graphRegularNeg = createTGraphErrors(negXValues, negPolarizations, negErrors, 20, 0.8, kBlack);

    setAxisLabelsAndRanges(graphRegularPos, "run index", yLabel, {0, static_cast<double>(xValues.size()) + 1}, yLimits);
    graphRegularPos->Draw("AP");
    graphRegularNeg->Draw("P SAME");

    // Draw the fitted constant lines on the right plot
    fitFuncPos->Draw("same");
    fitFuncNeg->Draw("same");

    // Draw outliers on the right plot
    if (graphPosOutliers) {
        graphPosOutliers->Draw("P SAME");
    }
    if (graphNegOutliers) {
        graphNegOutliers->Draw("P SAME");
    }

    // Create and draw text box in the top right for positive polarization
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.0275);
    text->DrawLatex(0.7, 0.7, Form("#mu_{+} = %.4g", muPos));
    text->DrawLatex(0.7, 0.65, Form("#sigma_{+} = %.4g", sigmaPos));
    text->DrawLatex(0.7, 0.6, Form("#chi^{2}/ndf_{+} = %.4g", chi2NdfPos));

    // Move text box to the left for negative polarization
    text->DrawLatex(0.225, 0.45, Form("#mu_{-} = %.4g", muNeg));
    text->DrawLatex(0.225, 0.4, Form("#sigma_{-} = %.4g", sigmaNeg));
    text->DrawLatex(0.225, 0.35, Form("#chi^{2}/ndf_{-} = %.4g", chi2NdfNeg));

    // Save the canvas to a file
    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs(outputFileName.c_str());

    // Clean up
    delete c;
    delete fitFuncPos;
    delete fitFuncNeg;
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

void addCanvasTopLabels(TCanvas* c, const std::vector<std::string>& q2_ranges) {
    c->cd(); // Switch to the entire canvas
    TLatex latex;
    latex.SetNDC();  // Use Normalized Device Coordinates
    latex.SetTextAlign(22);  // Centered alignment
    latex.SetTextSize(0.015);  // Adjust text size

    // Position and draw each Q² label
    for (size_t i = 0; i < q2_ranges.size(); ++i) {
        double xPos = 0.1 + i * 0.2;  // Adjust the 0.1 and 0.2 to space the labels evenly
        latex.DrawLatex(xPos, 0.99, q2_ranges[i].c_str());  // Position text above the plots
    }
}

void plotQ2yz_pT(const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData) {
    // Fit types, corresponding y-axis labels, and output file names
    std::vector<std::string> fitTypes = {"ALUsinphi", "AULoffset", "AULsinphi", "AULsin2phi", "ALL", "ALLcosphi", "doubleratio"};
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi}/F_{UU}",
        "A_{UL}^{offset}",
        "F_{UL}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin2#phi}/F_{UU}",
        "F_{LL}/F_{UU}",
        "F_{LL}^{cos#phi}/F_{UU}",
        "-F_{LU}^{sin#phi}/F_{LL}"
    };
    std::vector<std::string> outputFiles = {
        "output/epX_plots/Q2yz_pT_ALUsinphi.png",
        "output/epX_plots/Q2yz_pT_AULoffset.png",
        "output/epX_plots/Q2yz_pT_AULsinphi.png",
        "output/epX_plots/Q2yz_pT_AULsin2phi.png",
        "output/epX_plots/Q2yz_pT_ALL.png",
        "output/epX_plots/Q2yz_pT_ALLcosphi.png",
        "output/epX_plots/Q2yz_pT_doubleratio.png"
    };

    // Define different maxError thresholds for each fit type
    std::vector<double> maxErrors = {0.0275, 0.05, 0.0275, 0.05, 0.075, 0.05, 0.05}; // Add threshold for doubleratio

    // Define different y-axis ranges for each fit type
    std::vector<std::pair<double, double>> yRangesPerPlot = {
        {-0.09, 0.09},  // ALUsinphi
        {-0.199, 0.049},  // AULoffset
        {-0.099, 0.099},  // AULsinphi
        {-0.099, 0.099},  // AULsin2phi
        {-0.199, 0.599},  // ALL
        {-0.199, 0.199},  // ALLcosphi
        {-0.249, 0.249}     // doubleratio
    };

    // Define the legend once before the loop
    TLegend *legend = new TLegend(0.225, 0.225, 0.9, 0.9); // Adjust position and size of the legend box
    std::vector<std::string> zRanges = {
        "0.10 < z < 0.25",
        "0.25 < z < 0.35",
        "0.35 < z < 0.45",
        "0.45 < z < 0.55",
        "0.55 < z < 0.75"
    };
    std::vector<int> colors = {kBlack, kRed, kGreen, kBlue, kMagenta}; // Colors for z ranges

    // Add entries to the legend
    for (size_t zIndex = 0; zIndex < zRanges.size(); ++zIndex) {
        TGraph *dummyGraph = new TGraph();
        dummyGraph->SetMarkerColor(colors[zIndex]);
        dummyGraph->SetMarkerStyle(20); // Style of the marker
        dummyGraph->SetMarkerSize(1.5); // Size of the marker to make it more visible
        legend->AddEntry(dummyGraph, zRanges[zIndex].c_str(), "P");

        // Cast to TLegendEntry to set the text color
        TLegendEntry *entry = (TLegendEntry*)legend->GetListOfPrimitives()->Last();
        entry->SetTextColor(colors[zIndex]);
    }
    legend->SetTextSize(0.05); // Adjust text size if needed
    legend->SetFillColor(0);   // Make background transparent
    legend->SetLineColor(1);   // Add border

    // Loop over each fit type and generate the corresponding plot
    for (size_t fitIndex = 0; fitIndex < fitTypes.size(); ++fitIndex) {
        // Setup canvas for this fit type
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

        // Get the specific maxError and y-axis range for this fit type
        double maxError = maxErrors[fitIndex];
        std::pair<double, double> yRange = yRangesPerPlot[fitIndex];

        // Loop through each Q2 prefix and corresponding z prefixes
        for (size_t row = 0; row < Q2_prefixes.size(); ++row) {
            for (size_t q2Index = 0; q2Index < Q2_prefixes[row].size(); ++q2Index) {
                int padIndex = row * 5 + q2Index + 1;
                c->cd(padIndex); // Move to the appropriate pad in the canvas

                // Adjust margins
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

                // Loop through z prefixes
                for (size_t zIndex = 0; zIndex < z_prefixes.size(); ++zIndex) {
                    std::string key = Q2_prefixes[row][q2Index] + z_prefixes[zIndex] + "chi2Fits" + fitTypes[fitIndex];
                    auto it = asymmetryData.find(key);

                    if (it == asymmetryData.end()) {
                        continue;
                    }

                    const auto &data = it->second;

                    // Filter data by error
                    auto filteredData = filterDataByError(data, maxError);

                    // Skip this z-bin if all points are filtered out
                    if (filteredData.empty()) continue;

                    std::vector<double> x, y, yErr;
                    for (const auto &entry : filteredData) {
                        x.push_back(entry[0]);
                        y.push_back(entry[1]);
                        yErr.push_back(entry[2]);
                    }

                    TGraphErrors *graph = createTGraphErrors(x, y, yErr, 20, 0.8, colors[zIndex]);

                    // Set the custom y-axis range here
                    setCustomAxisLabelsAndRanges(graph, "P_{T} (GeV)", yLabels[fitIndex], {0.1, 0.9}, yRange);
                    
                    // Draw the graph
                    if (!firstGraphDrawn) {
                        graph->Draw("AP");
                        firstGraphDrawn = true;
                    } else {
                        graph->Draw("P SAME");
                    }

                    anyGraphDrawn = true;
                }

                if (!anyGraphDrawn) {
                    // Handle empty plot scenario
                    std::vector<double> dummyX = {-9999};
                    std::vector<double> dummyY = {0};
                    std::vector<double> dummyYErr = {0};
                    TGraphErrors *dummyGraph = createTGraphErrors(dummyX, dummyY, dummyYErr, 20, 0.8, kWhite);
                    setCustomAxisLabelsAndRanges(dummyGraph, "P_{T} (GeV)", yLabels[fitIndex], {0.1, 0.9}, yRange);
                    drawEmptyPlot(dummyGraph, q2Index, row, Q2_prefixes.size());
                }
                // Draw horizontal line except in certain positions
                if (!(row == 2 && q2Index == 4) && (row != 3 || (q2Index != 3 && q2Index != 4))) {
                    TLine *line = new TLine(0.15, 0.0, 0.95, 0.0);
                    line->SetLineColor(kGray + 2);
                    line->SetLineStyle(7);
                    line->Draw("same");
                }
            }
        }
        // Add y-range labels on the right-hand side
        addCanvasSideLabels(c, yRanges);

        // Add Q² labels at the top
        addCanvasTopLabels(c, topRowTitles);

        // Add legend to the canvas
        c->cd(20); // Navigate to the pad where the legend will be drawn
        legend->Draw();

        // Save the canvas to file
        gSystem->Exec("mkdir -p output/epX_plots");
        c->SaveAs(outputFiles[fitIndex].c_str());

        // Clean up canvas
        delete c;
    }

    // Clean up legend
    delete legend;
}

void plotQ2Dependence(
    const std::map<std::string, std::map<std::string, std::vector<std::vector<double>>>> &allVectors) {

    // Colors for different vector sets, avoiding black for better visibility
    std::vector<int> colors = {kRed, kBlue, kBlack, kGreen}; 
    std::vector<std::string> vectorNames = {"z1pT2y1", "z2pT2y1", "z1pT2y2", "z2pT2y2"};
    std::vector<std::string> legendLabels = {
        "0.65<y<0.75, 0.10<z<0.25, 0.21<PT<0.34", 
        "0.65<y<0.75, 0.25<z<0.35, 0.21<PT<0.34", 
        "0.55<y<0.65, 0.10<z<0.25, 0.21<PT<0.34", 
        "0.55<y<0.65, 0.25<z<0.35, 0.21<PT<0.34"
    };

    // Create a canvas with 1 row and 3 columns
    TCanvas *c = new TCanvas("c", "Q2 Dependence Plots", 2400, 800);
    c->Divide(3, 1);

    // Define the suffixes, y-axis labels, and ranges for each plot
    std::vector<std::string> suffixes = {"ALUsinphi", "ALL", "doubleratio"};
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi}/F_{UU}",
        "F_{LL}/F_{UU}",
        "-F_{LU}^{sin#phi}/F_{LL}"
    };
    std::vector<std::pair<double, double>> yRanges = {
        {-0.099, 0.019},  // ALUsinphi
        {-0.199, 0.599},  // ALL
        {-0.019, 0.099}   // doubleratio
    };
    
    std::string xLabel = "Q^{2} (GeV^{2})";
    std::pair<double, double> xLimits = {1.0, 7.0};

    // Loop through each suffix to create each plot
    for (size_t i = 0; i < suffixes.size(); ++i) {
        c->cd(i + 1); // Go to the i-th pad
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.15);

        // Adjusting legend position and font size
        double legendX1 = 0.18;  // Lower left corner X
        double legendY1 = 0.15;  // Lower left corner Y
        double legendX2 = 0.80;  // Upper right corner X
        double legendY2 = 0.35;  // Upper right corner Y
        TLegend *legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
        legend->SetTextSize(0.025); // Adjusted text size (decreased slightly)
        legend->SetBorderSize(1);

        // Loop through each vector set to plot them on the same canvas
        for (size_t j = 0; j < vectorNames.size(); ++j) {
            const std::string &vectorName = vectorNames[j];
            
            // Skip z2pT2y1 and z2pT2y2 for the doubleratio plot
            if (i == 2 && (vectorName == "z2pT2y1" || vectorName == "z2pT2y2")) {
                continue;
            }

            const auto &data = allVectors.at(vectorName).at(suffixes[i]);

            std::vector<double> x, y, yErr;
            for (const auto &point : data) {
                x.push_back(point[0]);     // Q² value
                y.push_back(point[1]);     // Asymmetry value
                yErr.push_back(point[2]);  // Asymmetry error
            }

            // Create a TGraphErrors and set axis labels and ranges
            TGraphErrors *graph = createTGraphErrors(x, y, yErr, 20, 0.8, colors[j]);
            setAxisLabelsAndRanges(graph, xLabel, yLabels[i], xLimits, yRanges[i]);

            // Draw the graph
            if (j == 0) {
                graph->Draw("AP");
            } else {
                graph->Draw("P SAME");
            }

            // Add the graph to the legend
            legend->AddEntry(graph, legendLabels[j].c_str(), "P");

            // Draw a horizontal line at y = 0
            TLine *line = new TLine(xLimits.first, 0, xLimits.second, 0);
            line->SetLineColor(kGray+2);
            line->SetLineStyle(7);
            line->Draw("same");
        }

        // Draw the legend on the current pad
        legend->Draw();
    }

    // Save the canvas to a file
    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs("output/epX_plots/Q2_dependence.png");

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
    plotDependence(asymmetryData, "x", "x_{B}", {0.06, 0.6}, "output/epX_plots/x_dependence_plots_comparison.png", "xall");
    // plotDependence(asymmetryData, "xall", "x_{B}", {0.06, 0.6}, "output/epX_plots/xall_dependence_plots.png");
    plotDependence(asymmetryData, "PT", "P_{T} (GeV)", {0.0, 1.0}, "output/epX_plots/PT_dependence_plots.png");
    plotDependence(asymmetryData, "PT", "P_{T} (GeV)", {0.0, 1.0}, "output/epX_plots/PT_dependence_plots_comparison.png", "PTall");
    // plotDependence(asymmetryData, "PTall", "P_{T} (GeV)", {0.0, 1.0}, "output/epX_plots/PTall_dependence_plots.png");
    plotDependence(asymmetryData, "xF", "x_{F}", {-0.8, 0.6}, "output/epX_plots/xF_dependence_plots.png");
    plotDependence(asymmetryData, "xF", "x_{F}", {-0.8, 0.6}, "output/epX_plots/xF_dependence_plots_comparison.png", "xFall");
    // plotDependence(asymmetryData, "xFall", "x_{F}", {-0.8, 0.6}, "output/epX_plots/xFall_dependence_plots.png");
    plotDependence(asymmetryData, "Mx", "M_{x} (GeV)", {0, 3}, "output/epX_plots/Mx_dependence_plots.png");
    // plotDependence(asymmetryData, "runnum", "run number", {16135, 16774}, "output/epX_plots/runnum_dependence_plots.png");
    // plotRunnumDependence(asymmetryData, "runnum", "run number", "output/epX_plots/runnum_dependence_plots.png");

    // // Plot PT and xF dependence comparison
    // plotComparison(asymmetryData, "output/epX_plots/PT_xF_dependence_comparison.png");
    // // Plot Q2-y-z dependence
    // plotQ2yz_pT(asymmetryData);

    // // Extract Q² dependence vectors
    // auto allVectors = extractQ2Dependence(asymmetryData, kinematicData);
    // plotQ2Dependence(allVectors);

    plotTargetPolarizationDependence(targetPolarizationData, "run number", "output/epX_plots/target_polarization_dependence.png");

    return 0;
}
