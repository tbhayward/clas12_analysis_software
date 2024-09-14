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

// good runs
std::vector<std::tuple<int, double, double>> targetPolarizationData = {{{16137, 0.7677, 0.0184}, {16138, 0.7226, 0.0179}, {16144, 0.7484, 0.0217}, {16145, 0.789, 0.0223}, {16146, 0.7063, 0.025}, {16148, 0.7459, 0.0219}, {16156, -0.6508, 0.0519}, {16157, -0.7146, 0.0293}, {16158, -0.7163, 0.0321}, {16164, -0.645, 0.1282}, {16166, -0.779, 0.0375}, {16167, -0.7985, 0.0265}, {16168, -0.779, 0.0788}, {16169, -0.7351, 0.0531}, {16170, -0.674, 0.0211}, {16178, -0.7082, 0.0247}, {16211, 0.8711, 0.0239}, {16214, 0.8096, 0.0231}, {16221, 0.8168, 0.023}, {16222, 0.808, 0.0217}, {16223, 0.7966, 0.0209}, {16224, 0.7067, 0.0682}, {16225, 0.7696, 0.0327}, {16226, 0.893, 0.0365}, {16228, 0.8948, 0.0284}, {16231, -0.7252, 0.0271}, {16232, -0.731, 0.0209}, {16233, -0.7518, 0.0205}, {16235, -0.6254, 0.0528}, {16236, -0.6591, 0.0337}, {16238, -0.5733, 0.0488}, {16243, -0.7715, 0.0244}, {16244, -0.7191, 0.0267}, {16245, -0.7029, 0.0196}, {16246, -0.6532, 0.0195}, {16248, -0.6896, 0.0205}, {16249, -0.6618, 0.0394}, {16250, -0.6346, 0.0446}, {16251, -0.8599, 0.0379}, {16252, -0.7241, 0.0331}, {16253, -0.8366, 0.0242}, {16256, -0.6868, 0.0203}, {16257, -0.696, 0.0202}, {16259, -0.6676, 0.0206}, {16260, -0.7266, 0.0254}, {16318, 0.7732, 0.023}, {16320, 0.7917, 0.0216}, {16321, 0.749, 0.0223}, {16322, 0.7914, 0.0285}, {16323, 0.8742, 0.0365}, {16325, 0.677, 0.0643}, {16326, 0.8518, 0.0455}, {16327, 0.7927, 0.0209}, {16328, 0.766, 0.0213}, {16329, 0.7746, 0.0204}, {16330, 0.78, 0.02}, {16331, 0.7265, 0.0293}, {16332, 0.7157, 0.0687}, {16333, 0.7987, 0.0435}, {16335, -0.7002, 0.0214}, {16336, -0.7016, 0.0266}, {16337, -0.6994, 0.0295}, {16338, -0.7127, 0.0222}, {16339, -0.8082, 0.0361}, {16341, -0.6767, 0.0208}, {16343, -0.7024, 0.0213}, {16345, -0.5727, 0.0329}, {16346, -0.7184, 0.0203}, {16348, -0.7732, 0.0211}, {16350, -0.7161, 0.0206}, {16352, -0.7331, 0.0212}, {16353, -0.6895, 0.0199}, {16354, -0.7316, 0.0209}, {16355, -0.7127, 0.0207}, {16356, -0.6927, 0.0208}, {16357, -0.6967, 0.024}, {16709, 0.6617, 0.0334}, {16710, 0.8414, 0.0314}, {16711, 0.7551, 0.0205}, {16712, 0.7663, 0.0204}, {16713, 0.7704, 0.0204}, {16715, 0.7103, 0.0247}, {16716, 0.7885, 0.0268}, {16717, 0.7316, 0.0216}, {16718, 0.6173, 0.0345}, {16719, 0.6784, 0.0205}, {16720, 0.7184, 0.0251}, {16721, -0.4516, 0.068}, {16722, -0.7295, 0.048}, {16723, -0.7011, 0.0219}, {16726, -0.6371, 0.0188}, {16727, -0.6975, 0.0211}, {16728, -0.6279, 0.0195}, {16729, -0.6943, 0.0388}, {16730, -0.6534, 0.0236}, {16731, -0.7173, 0.0276}, {16732, -0.6901, 0.0244}, {16733, -0.7308, 0.0484}, {16734, -0.5724, 0.0968}, {16736, -0.6599, 0.02}, {16738, -0.6721, 0.0239}, {16743, -0.7407, 0.0215}, {16744, -0.6661, 0.0191}, {16746, -0.6751, 0.0196}, {16747, -0.6623, 0.0193}, {16748, -0.6544, 0.0198}, {16749, -0.6525, 0.0438}, {16750, -0.7157, 0.0201}, {16751, -0.6532, 0.0378}, {16752, -0.6773, 0.0301}, {16754, -0.6698, 0.0194}, {16755, -0.6767, 0.0195}, {16756, -0.6443, 0.0377}, {16757, -0.6146, 0.0191}, {16759, -0.6446, 0.0187}, {16761, -0.6492, 0.0206}, {16762, -0.6736, 0.0198}, {16763, -0.6451, 0.0244}, {16765, -0.6837, 0.0212}, {16766, -0.6931, 0.0298}, {16767, 0.8475, 0.0223}, {16768, 0.7774, 0.0212}, {16769, 0.7328, 0.0272}, {16770, 0.7803, 0.0219}, {16771, 0.7577, 0.0204}, {16772, 0.7756, 0.022}}};
// good and bad runs
// std::vector<std::tuple<int, double, double>> targetPolarizationData = {
//     {16137, 0.7677, 0.0184}, {16138, 0.7226, 0.0179}, {16144, 0.7484, 0.0217},
//     {16145, 0.789, 0.0223}, {16146, 0.7063, 0.025}, {16148, 0.7459, 0.0219},
//     {16156, -0.6508, 0.0519}, {16157, -0.7146, 0.0293}, {16158, -0.7163, 0.0321},
//     {16164, -0.645, 0.1282}, {16166, -0.779, 0.0375}, {16167, -0.7985, 0.0265},
//     {16168, -0.779, 0.0788}, {16169, -0.7351, 0.0531}, {16170, -0.674, 0.0211},
//     {16178, -0.7082, 0.0247}, {16211, 0.8711, 0.0239}, {16213, 0.92, 0.028},
//     {16214, 0.8096, 0.0231}, {16221, 0.8168, 0.023}, {16222, 0.808, 0.0217},
//     {16223, 0.7966, 0.0209}, {16224, 0.7067, 0.0682}, {16225, 0.7696, 0.0327},
//     {16226, 0.893, 0.0365}, {16228, 0.8948, 0.0284}, {16231, -0.7252, 0.0271},
//     {16232, -0.731, 0.0209}, {16233, -0.7518, 0.0205}, {16234, -0.9052, 0.0398},
//     {16235, -0.6254, 0.0528}, {16236, -0.6591, 0.0337}, {16238, -0.5733, 0.0488},
//     {16243, -0.7715, 0.0244}, {16244, -0.7191, 0.0267}, {16245, -0.7029, 0.0196},
//     {16246, -0.6532, 0.0195}, {16248, -0.6896, 0.0205}, {16249, -0.6618, 0.0394},
//     {16250, -0.6346, 0.0446}, {16251, -0.8599, 0.0379}, {16252, -0.7241, 0.0331},
//     {16253, -0.8366, 0.0242}, {16256, -0.6868, 0.0203}, {16257, -0.696, 0.0202},
//     {16259, -0.6676, 0.0206}, {16260, -0.7266, 0.0254}, {16317, 0.8632, 0.0433},
//     {16318, 0.7732, 0.023}, {16320, 0.7917, 0.0216}, {16321, 0.749, 0.0223},
//     {16322, 0.7914, 0.0285}, {16323, 0.8742, 0.0365}, {16325, 0.677, 0.0643},
//     {16326, 0.8518, 0.0455}, {16327, 0.7927, 0.0209}, {16328, 0.766, 0.0213},
//     {16329, 0.7746, 0.0204}, {16330, 0.78, 0.02}, {16331, 0.7265, 0.0293},
//     {16332, 0.7157, 0.0687}, {16333, 0.7987, 0.0435}, {16335, -0.7002, 0.0214},
//     {16336, -0.7016, 0.0266}, {16337, -0.6994, 0.0295}, {16338, -0.7127, 0.0222},
//     {16339, -0.8082, 0.0361}, {16341, -0.6767, 0.0208}, {16343, -0.7024, 0.0213},
//     {16345, -0.5727, 0.0329}, {16346, -0.7184, 0.0203}, {16348, -0.7732, 0.0211},
//     {16350, -0.7161, 0.0206}, {16352, -0.7331, 0.0212}, {16353, -0.6895, 0.0199},
//     {16354, -0.7316, 0.0209}, {16355, -0.7127, 0.0207}, {16356, -0.6927, 0.0208},
//     {16357, -0.6967, 0.024}, {16658, 0.3326, 0.0276}, {16659, 0.4347, 0.0257},
//     {16660, 0.642, 0.0311}, {16664, 0.5856, 0.0562}, {16665, 0.5405, 0.0184},
//     {16666, 0.5042, 0.0249}, {16671, 0.6526, 0.0205}, {16672, 0.5568, 0.0198},
//     {16673, 0.5191, 0.0473}, {16674, 0.629, 0.0326}, {16675, 0.5129, 0.0194},
//     {16676, 0.5079, 0.0182}, {16678, 0.5775, 0.0284}, {16679, 0.5236, 0.0271},
//     {16681, 0.6906, 0.0195}, {16682, 0.7089, 0.0202}, {16683, 0.578, 0.0188},
//     {16685, 0.5514, 0.0393}, {16686, 0.6273, 0.0238}, {16687, 0.5774, 0.1147},
//     {16688, 0.6294, 0.0189}, {16689, 0.8544, 0.0351}, {16690, 0.6241, 0.0491},
//     {16692, 0.6865, 0.0211}, {16693, 0.7503, 0.0204}, {16695, 0.2715, 0.0911},
//     {16709, 0.6617, 0.0334}, {16710, 0.8414, 0.0314}, {16711, 0.7551, 0.0205},
//     {16712, 0.7663, 0.0204}, {16713, 0.7704, 0.0204}, {16715, 0.7103, 0.0247},
//     {16716, 0.7885, 0.0268}, {16717, 0.7316, 0.0216}, {16718, 0.6173, 0.0345},
//     {16719, 0.6784, 0.0205}, {16720, 0.7184, 0.0251}, {16721, -0.4516, 0.068},
//     {16722, -0.7295, 0.048}, {16723, -0.7011, 0.0219}, {16726, -0.6371, 0.0188},
//     {16727, -0.6975, 0.0211}, {16728, -0.6279, 0.0195}, {16729, -0.6943, 0.0388},
//     {16730, -0.6534, 0.0236}, {16731, -0.7173, 0.0276}, {16732, -0.6901, 0.0244},
//     {16733, -0.7308, 0.0484}, {16734, -0.5724, 0.0968}, {16736, -0.6599, 0.02},
//     {16738, -0.6721, 0.0239}, {16742, -0.655, 0.0354}, {16743, -0.7407, 0.0215},
//     {16744, -0.6661, 0.0191}, {16746, -0.6751, 0.0196}, {16747, -0.6623, 0.0193},
//     {16748, -0.6544, 0.0198}, {16749, -0.6525, 0.0438}, {16750, -0.7157, 0.0201},
//     {16751, -0.6532, 0.0378}, {16752, -0.6773, 0.0301}, {16753, -0.8703, 0.0337},
//     {16754, -0.6698, 0.0194}, {16755, -0.6767, 0.0195}, {16756, -0.6443, 0.0377},
//     {16757, -0.6146, 0.0191}, {16759, -0.6446, 0.0187}, {16761, -0.6492, 0.0206},
//     {16762, -0.6736, 0.0198}, {16763, -0.6451, 0.0244}, {16765, -0.6837, 0.0212},
//     {16766, -0.6931, 0.0298}, {16767, 0.8475, 0.0223}, {16768, 0.7774, 0.0212},
//     {16769, 0.7328, 0.0272}, {16770, 0.7803, 0.0219}, {16771, 0.7577, 0.0204},
//     {16772, 0.7756, 0.022}
// };

// Global variables for systematic uncertainties
const double LU_SYS_UNCERTAINTY = 0.029;
const double UL_SYS_UNCERTAINTY = 0.101;
const double LL_SYS_UNCERTAINTY_LU_COMPONENT = LU_SYS_UNCERTAINTY;
const double LL_SYS_UNCERTAINTY_UL_COMPONENT = UL_SYS_UNCERTAINTY;



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
                    ratioValue=ratioValue;
                    doubleratioData.push_back({xValue, ratioValue, error});
                }

                // Store the doubleratio data
                asymmetryData[baseKey + "doubleratio"] = doubleratioData;
            }
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


void plotDependence(
    const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData,
    const std::string &prefix, 
    const std::string &xLabel, 
    const std::pair<double, double> &xLimits, 
    const std::string &outputFileName,
    const std::string &extraPrefix = "") {  // Optional extra prefix
    TCanvas *c = new TCanvas("c", "Dependence Plots", 1200, 800);
    c->Divide(3, 2);

    // Updated list of suffixes (removing AULoffset and adding doubleratio)
    std::vector<std::string> suffixes = {"ALUsinphi", "AULsinphi", "AULsin2phi", "ALL", "doubleratio", "ALLcosphi"};
    
    // Corresponding y-axis labels
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin(2#phi)}/F_{UU}",
        "F_{LL}/F_{UU}",
        "-F_{LU}^{sin#phi}/F_{LL}", 
        "F_{LL}^{cos#phi}/F_{UU}"
    };

    // Updated loop to accommodate the new plot order
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
                                   (suffixes[i] == "ALL") ? std::make_pair(0.15, 0.4) :
                                   (suffixes[i] == "doubleratio") ? std::make_pair(-0.1, 0.2) :
                                   std::make_pair(0.00, 0.05));  // Adjusted y-axis range for doubleratio
            graphStat->Draw("AP");

            // Add the dashed gray line at y = 0
            TLine *line = new TLine(xLimits.first, 0, xLimits.second, 0);
            line->SetLineColor(kGray+2);
            line->SetLineStyle(7);  // Dashed line
            line->Draw("same");

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

                if (i == 3 || i == 4) {  
                    // Add the text box for labels only if the extraPrefix is provided
                    TPaveText *text = new TPaveText(0.18, 0.7, 0.43, 0.9, "NDC");
                    text->SetTextAlign(13);
                    text->SetBorderSize(1);  // Set border size to 1 for a black border
                    text->SetFillColor(0);
                    text->AddText("#font[42]{M_{x} > 1.35 GeV}");  // Black text line
                    text->AddText("#font[42]{#color[2]{M_{x} > 0.55 GeV}}");  // Red text line
                    text->Draw();
                } else {
                    // Add the text box for labels only if the extraPrefix is provided
                    TPaveText *text = new TPaveText(0.65, 0.7, 0.9, 0.9, "NDC");
                    text->SetTextAlign(13);
                    text->SetBorderSize(1);  // Set border size to 1 for a black border
                    text->SetFillColor(0);
                    text->AddText("#font[42]{M_{x} > 1.35 GeV}");  // Black text line
                    text->AddText("#font[42]{#color[2]{M_{x} > 0.55 GeV}}");  // Red text line
                    text->Draw();
                }
            }
        }
    }

    gSystem->Exec("mkdir -p output/rho0_plots");
    c->SaveAs(outputFileName.c_str());
    delete c;
}

void plotCombinationDependence(
    const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData,
    const std::string &prefix0,  // Baseline dataset (gray)
    const std::string &prefix1,  // Dataset 1 (orange)
    const std::string &prefix2,  // Dataset 2 (blue)
    const std::string &xLabel, 
    const std::pair<double, double> &xLimits, 
    const std::pair<double, double> &yRangeALU,  // y range for ALUsinphi
    const std::pair<double, double> &yRangeALL,  // y range for ALL
    const std::string &outputFileName,
    const std::vector<std::string> &legendEntries) {  // Should contain 3 entries now

    // Create the canvas and divide it into 2 subplots (1x2)
    TCanvas *c = new TCanvas("c", "Combination Dependence Plots", 1800, 600); // Adjust the canvas size for a 1x2 layout
    c->Divide(2, 1);  // 1 row and 2 columns

    // Loop over the two subplots
    for (size_t i = 0; i < 2; ++i) {
        c->cd(i + 1);  // Move to the appropriate pad

        if (i == 0) {
            // Left subplot: reduce the right margin to minimize space between plots
            gPad->SetLeftMargin(0.05);
            gPad->SetRightMargin(0.05);  // Decrease right margin
            gPad->SetBottomMargin(0.05);
            gPad->SetTopMargin(0.05);
        } else {
            // Right subplot: reduce the left margin to minimize space between plots
            gPad->SetLeftMargin(0.05);  // Decrease left margin
            gPad->SetRightMargin(0.05);  // Standard right margin
            gPad->SetBottomMargin(0.05);
            gPad->SetTopMargin(0.05);
        }
    }

    // Define the two suffixes and corresponding y-axis labels
    std::vector<std::string> suffixes = {"ALUsinphi", "ALL"};
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi_{#pi^{+}}}/F_{UU}",
        "F_{LL}/F_{UU}"
    };

    // Declare graph0, graph1, and graph2 outside of the loop so that they can be used for the legend
    TGraphErrors *graph0 = nullptr;
    TGraphErrors *graph1 = nullptr;
    TGraphErrors *graph2 = nullptr;

    // Loop over the two suffixes to create the subplots
    for (size_t i = 0; i < suffixes.size(); ++i) {
        c->cd(i + 1);  // Move to the appropriate pad
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.15);

        // Build the keys for all three datasets
        std::string key0 = prefix0 + "chi2Fits" + suffixes[i];
        std::string key1 = prefix1 + "chi2Fits" + suffixes[i];
        std::string key2 = prefix2 + "chi2Fits" + suffixes[i];

        // Check if all three datasets exist
        auto it0 = asymmetryData.find(key0);
        auto it1 = asymmetryData.find(key1);
        auto it2 = asymmetryData.find(key2);

        if (it0 != asymmetryData.end() && it1 != asymmetryData.end() && it2 != asymmetryData.end()) {
            const auto &data0 = it0->second;
            const auto &data1 = it1->second;
            const auto &data2 = it2->second;

            // Extract values for the baseline dataset (prefix0)
            std::vector<double> x0, y0, y0Err;
            for (const auto &entry : data0) {
                x0.push_back(entry[0]);
                y0.push_back(entry[1]);
                y0Err.push_back(entry[2]);
            }

            // Extract values for the first dataset (prefix1)
            std::vector<double> x1, y1, y1Err;
            for (const auto &entry : data1) {
                x1.push_back(entry[0]);
                y1.push_back(entry[1]);
                y1Err.push_back(entry[2]);
            }

            // Extract values for the second dataset (prefix2)
            std::vector<double> x2, y2, y2Err;
            for (const auto &entry : data2) {
                x2.push_back(entry[0]);
                y2.push_back(entry[1]);
                y2Err.push_back(entry[2]);
            }

            // Create TGraphErrors for all three datasets (using gray, orange, and blue circles)
            graph0 = createTGraphErrors(x0, y0, y0Err, 20, 0.8, kGray+5);   // Gray circles
            graph1 = createTGraphErrors(x1, y1, y1Err, 20, 0.8, kOrange+5); // Orange circles
            graph2 = createTGraphErrors(x2, y2, y2Err, 20, 0.8, kBlue);   // Blue circles

            // Set axis labels and ranges for the graph
            setAxisLabelsAndRanges(graph0, xLabel, yLabels[i], xLimits, 
                                   (suffixes[i] == "ALL") ? yRangeALL : yRangeALU);

            // Draw all three graphs on the same pad
            graph0->Draw("AP");  // Baseline (gray)
            graph1->Draw("P SAME");  // First dataset (orange)
            graph2->Draw("P SAME");  // Second dataset (blue)

            // Add a dashed gray line at y = 0
            TLine *line = new TLine(xLimits.first, 0, xLimits.second, 0);
            line->SetLineColor(kGray+2);
            line->SetLineStyle(7);  // Dashed line
            line->Draw("same");

            // Create a legend for each subplot, positioned at the top right
            TLegend *legend = new TLegend(0.35, 0.75, 0.95, 0.95);  // Set the fixed position you wanted
            legend->SetBorderSize(1);  // Set border size to 1 for a black border
            legend->SetTextSize(0.0325);  // Set smaller text size

            // Entry for prefix0 (gray)
            legend->AddEntry(graph0, legendEntries[0].c_str(), "p");
            TLegendEntry *entry0 = dynamic_cast<TLegendEntry*>(legend->GetListOfPrimitives()->Last());  // Get last entry and cast it to TLegendEntry
            if (entry0) entry0->SetTextColor(kGray);  // Set color of the first entry to gray

            // Entry for prefix1 (orange)
            legend->AddEntry(graph1, legendEntries[1].c_str(), "p");
            TLegendEntry *entry1 = dynamic_cast<TLegendEntry*>(legend->GetListOfPrimitives()->Last());  // Get last entry and cast it to TLegendEntry
            if (entry1) entry1->SetTextColor(kOrange);  // Set color of the second entry to orange

            // Entry for prefix2 (blue)
            legend->AddEntry(graph2, legendEntries[2].c_str(), "p");
            TLegendEntry *entry2 = dynamic_cast<TLegendEntry*>(legend->GetListOfPrimitives()->Last());  // Get last entry and cast it to TLegendEntry
            if (entry2) entry2->SetTextColor(kBlue);  // Set color of the third entry to blue

            legend->Draw();
        }
    }

    // Save the canvas as a PNG file
    gSystem->Exec("mkdir -p output/rho0_plots");
    c->SaveAs(outputFileName.c_str());

    // Clean up memory
    delete c;
}


int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <asymmetries.txt>\n";
        return 1;
    }

    std::string asymmetryFile = argv[1];

    // Read the asymmetry data from the file
    std::map<std::string, std::vector<std::vector<double>>> asymmetryData = readAsymmetries(asymmetryFile);

    // // Print out the parsed data
    // std::cout << "Asymmetry Data:\n";
    // printData(asymmetryData);

    // Call the plotting function for different dependencies
    // plotDependence(asymmetryData, "epiplus", "P_{T} (GeV)", {0.0, 1.1}, "output/rho0_plots/PT_epiplus_dependence_plots.png");
    // plotDependence(asymmetryData, "epipluspiminus", "P_{T} (GeV)", {0.0, 1.1}, "output/rho0_plots/PT_epipluspiminus_dependence_plots.png");
    // plotDependence(asymmetryData, "epipluspiminus_rho0_free", "P_{T} (GeV)", {0.0, 1.1}, "output/rho0_plots/PT_epipluspiminus_rho0_free_dependence_plots.png");
    // plotDependence(asymmetryData, "eppiplus", "P_{T} (GeV)", {0.0, 1.1}, "output/rho0_plots/PT_eppiplus_dependence_plots.png");
    // plotDependence(asymmetryData, "eppiplus_rho0_free", "P_{T} (GeV)", {0.0, 1.1}, "output/rho0_plots/PT_eppiplus_rho0_free_dependence_plots.png");
    // plotDependence(asymmetryData, "eppipluspiminus", "P_{T} (GeV)", {0.0, 1.1}, "output/rho0_plots/PT_eppipluspiminus_dependence_plots.png");
    // plotDependence(asymmetryData, "eppipluspiminus_rho0_free_A", "P_{T} (GeV)", {0.0, 1.1}, "output/rho0_plots/PT_eppipluspiminus_rho0_free_A_dependence_plots.png");
    // plotDependence(asymmetryData, "eppipluspiminus_rho0_free_B", "P_{T} (GeV)", {0.0, 1.1}, "output/rho0_plots/PT_eppipluspiminus_rho0_free_B_dependence_plots.png");
  
    // Plot combination for epiplus, epipluspiminus, and epipluspiminus_rho0_free
    plotCombinationDependence(asymmetryData, 
        "epiplus",  // Baseline dataset (black)
        "epipluspiminus",  // First dataset (red)
        "epipluspiminus_rho0_free",  // Second dataset (blue)
        "P_{T} (GeV)", 
        {0.0, 1.1}, 
        {0.0, 0.035},  // y range for ALUsinphi
        {0.2, 0.4},  // y range for ALL
        "output/rho0_plots/PT_epipluspiminus_combination_dependence_plots.png", 
        {
            "e'#pi^{+}X, M_{x (#pi^{+})} > 1.5 (GeV)",  // Baseline (black)
            "e'#pi^{+}#pi^{-}X, M_{x (#pi^{+})} > 1.5 (GeV)",  // Dataset 1 (red)
            "e'#pi^{+}#pi^{-}X, M_{x (#pi^{+})} > 1.5 (GeV), M_{x (#pi^{+}#pi^{-})} > 1.05 (GeV)"  // Dataset 2 (blue)
        }
    );

    return 0;
}
