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

    // Updated list of suffixes (removing AULoffset and adding doubleratio)
    std::vector<std::string> suffixes = {"ALUsinphi", "AULsinphi", "AULsin2phi", "ALL", "doubleratio", "ALLcosphi"};
    
    // Corresponding y-axis labels
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin(2#phi)}/F_{UU}",
        "F_{LL}/F_{UU}",
        "-F_{LU}^{sin#phi}/F_{LL}",  // yLabel for doubleratio
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
                                   (suffixes[i] == "ALL") ? std::make_pair(-0.1, 0.6) :
                                   (suffixes[i] == "doubleratio") ? std::make_pair(-0.02, 0.3) :
                                   std::make_pair(-0.1, 0.1));  // Adjusted y-axis range for doubleratio
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

    // Create a canvas with 1 row and 2 columns for original plots
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
            if (std::abs(polarization - muPos) > 3.5 * error) {
                outlierPosX.push_back(xValues[i]);
                outlierPosY.push_back(polarization);
                outlierPosErrors.push_back(error);
                std::cout << "Outlier (Positive) found: Run Number " << runNumbers[i] << ", Deviation: " << std::abs(polarization - muPos) / error << std::endl;
            }
        } else {
            if (std::abs(polarization - muNeg) > 3.5 * error) {
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

    // Save the first canvas
    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs(outputFileName.c_str());

    // Create a single canvas and plot both histograms together
    TCanvas *c2 = new TCanvas("c2", "Target Polarization Histograms", 1600, 800);

    // Create histograms for positive and negative target polarizations
    TH1F *histPos = new TH1F("histPos", "Positive Target Polarizations", 50, 0.0, 1.0);  // 0 to 1 for positives
    TH1F *histNeg = new TH1F("histNeg", "Negative Target Polarizations", 50, -1.0, 0.0); // -1 to 0 for negatives

    // Remove stat boxes
    gStyle->SetOptStat(0);

    // Fill histograms with positive and negative polarization values
    for (double val : posPolarizations) {
        histPos->Fill(val);
    }
    for (double val : negPolarizations) {
        histNeg->Fill(val);
    }

    // Fit the histograms to a Gaussian distribution
    TF1 *gausPos = new TF1("gausPos", "gaus", 0.0, 1.0);  // 0 to 1 for positives
    TF1 *gausNeg = new TF1("gausNeg", "gaus", -1.0, 0.0); // -1 to 0 for negatives

    histPos->Fit(gausPos, "Q");  // Silent mode for positive fit
    histNeg->Fit(gausNeg, "Q");  // Silent mode for negative fit

    // Retrieve fit parameters for the legend
    double muPosHist = gausPos->GetParameter(1);
    double sigmaPosHist = gausPos->GetParameter(2);
    double muNegHist = gausNeg->GetParameter(1);
    double sigmaNegHist = gausNeg->GetParameter(2);

    // Set histogram line styles
    histPos->SetLineColor(kRed);
    histNeg->SetLineColor(kBlue);

    // Set the x-axis range to cover both positive and negative target polarizations
    TH1F *frame = new TH1F("frame", "", 100, -1.0, 1.0);  // Frame to set the range of the plot
    frame->SetMinimum(0); // Minimum y value to ensure the histograms are fully visible
    frame->SetMaximum(25); // Adjust this as per your data if necessary
    frame->GetXaxis()->SetTitle("Target Polarization");
    frame->GetYaxis()->SetTitle("Runs");
    frame->Draw();  // Draw the frame first

    // Draw positive polarization histogram as red data points without horizontal error bars
    histPos->SetMarkerColor(kRed);
    histPos->SetMarkerStyle(20);
    histPos->SetMarkerSize(1.0);
    histPos->Draw("E SAME");  // Draw histogram with vertical error bars only

    // Draw fitted Gaussian curve for positive histogram
    gausPos->SetLineColor(kRed);
    gausPos->Draw("SAME");

    // Draw negative polarization histogram as blue data points without horizontal error bars
    histNeg->SetMarkerColor(kBlue);
    histNeg->SetMarkerStyle(21);
    histNeg->SetMarkerSize(1.0);
    histNeg->Draw("E SAME");  // Draw histogram with vertical error bars only

    // Draw fitted Gaussian curve for negative histogram
    gausNeg->SetLineColor(kBlue);
    gausNeg->Draw("SAME");

    // Add a legend for the fits and corresponding means and sigmas
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(gausPos, Form("#mu_{+} = %.4f, #sigma_{+} = %.4f ", muPosHist, sigmaPosHist), "l")->SetTextColor(kRed);
    leg->AddEntry(gausNeg, Form("#mu_{-} = %.4f, #sigma_{-} = %.4f ", muNegHist, sigmaNegHist), "l")->SetTextColor(kBlue);
    leg->Draw();

    // Save the canvas
    c2->SaveAs("output/epX_plots/target_polarization_histograms_combined.png");

    // Clean up
    delete c2;
    delete histPos;
    delete histNeg;
    delete gausPos;
    delete gausNeg;
    delete leg;
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

void plotMultipleQ2multiDependence(
    const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData,
    const std::vector<std::string> &prefixes,  // List of prefixes
    const std::string &xLabel, 
    const std::pair<double, double> &xLimits, 
    const std::string &outputFileName) {

    TCanvas *c = new TCanvas("c", "Q2multi Dependence Plots", 1200, 800);
    c->Divide(3, 2);

    // Updated list of suffixes (removing AULoffset and adding doubleratio)
    std::vector<std::string> suffixes = {"ALUsinphi", "AULsinphi", "AULsin2phi", "ALL", "doubleratio", "ALLcosphi"};
    
    // Corresponding y-axis labels
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin(2#phi)}/F_{UU}",
        "F_{LL}/F_{UU}",
        "-F_{LU}^{sin#phi}/F_{LL}",  // yLabel for doubleratio
        "F_{LL}^{cos#phi}/F_{UU}"
    };

    // Colors and marker styles for each prefix
    std::vector<int> colors = {kRed, kBlue, kGreen};  // Updated colors
    std::vector<int> markers = {20, 21, 22};  // Circle, Square, Triangle markers

    // Text labels for the legend
    std::vector<std::string> legendLabels = {
        "(0.12 < x < 0.15)", 
        "(0.15 < x < 0.18)", 
        "(0.18 < x < 0.21)"
    };

    // Title to add to each subplot
    std::string plotTitle = "0.16 < z < 0.22, 0.325 < P_{T} < 0.475, M_{x} > 0.95";

    // Updated loop to accommodate the new plot order
    for (size_t i = 0; i < suffixes.size(); ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.15);

        bool firstGraphDrawn = false;

        // Create a legend in the top right with a border and background
        TLegend *legend = new TLegend(0.55, 0.7, 0.9, 0.9);  // Adjust position and size
        legend->SetTextSize(0.035);  // Adjust text size
        legend->SetBorderSize(1);  // Set border size
        legend->SetFillStyle(1001);   // Solid white background

        for (size_t p = 0; p < prefixes.size(); ++p) {
            std::string key = prefixes[p] + "chi2Fits" + suffixes[i];
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

                TGraphErrors *graphStat = createTGraphErrors(x, y, yStatErr, markers[p], 1.0, colors[p]);  // Adjusted marker size to 1.0
                graphStat->SetLineWidth(1);  // Ensure the line is centered properly
                
                setAxisLabelsAndRanges(graphStat, xLabel, yLabels[i], xLimits, 
                                       (suffixes[i] == "ALL") ? std::make_pair(-0.1, 0.5) :
                                       (suffixes[i] == "doubleratio") ? std::make_pair(-0.02, 0.4) :
                                       std::make_pair(-0.06, 0.06));  // Adjusted y-axis range for doubleratio
                
                if (!firstGraphDrawn) {
                    graphStat->Draw("AP");
                    firstGraphDrawn = true;
                } else {
                    graphStat->Draw("P SAME");
                }

                // Add each entry to the legend with the corresponding color and label
                TLegendEntry* legendEntry = legend->AddEntry(graphStat, legendLabels[p].c_str(), "p");
                legendEntry->SetTextColor(colors[p]);  // Set the text color to match the graph color
            }
        }

        // Add the dashed gray line at y = 0
        TLine *line = new TLine(xLimits.first, 0, xLimits.second, 0);
        line->SetLineColor(kGray+2);
        line->SetLineStyle(7);  // Dashed line
        line->Draw("same");

        // Draw the legend in the current subplot
        legend->Draw();

        // Add a title to each subplot
        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.035);
        latex->SetTextAlign(22);  // Centered text alignment
        latex->DrawLatex(0.535, 0.93, plotTitle.c_str());  // Adjust position and size
    }

    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs(outputFileName.c_str());
    delete c;
}

void plotMultipleQ2extramultiDependence(
    const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData,
    const std::vector<std::string> &prefixes,  // List of prefixes
    const std::string &xLabel, 
    const std::pair<double, double> &xLimits, 
    const std::string &outputFileName) {

    TCanvas *c = new TCanvas("c", "Q2extramulti Dependence Plots", 1200, 800);
    c->Divide(3, 2);

    // Updated list of suffixes (removing AULoffset and adding doubleratio)
    std::vector<std::string> suffixes = {"ALUsinphi", "AULsinphi", "AULsin2phi", "ALL", "doubleratio", "ALLcosphi"};
    
    // Corresponding y-axis labels
    std::vector<std::string> yLabels = {
        "F_{LU}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin#phi}/F_{UU}",
        "F_{UL}^{sin(2#phi)}/F_{UU}",
        "F_{LL}/F_{UU}",
        "-F_{LU}^{sin#phi}/F_{LL}",  // yLabel for doubleratio
        "F_{LL}^{cos#phi}/F_{UU}"
    };

    // Colors and marker styles for each prefix
    std::vector<int> colors = {kRed, kBlue, kGreen};  // Updated colors
    std::vector<int> markers = {20, 21, 22};  // Circle, Square, Triangle markers

    // Text labels for the legend
    std::vector<std::string> legendLabels = {
        "(0.12 < x < 0.15)", 
        "(0.15 < x < 0.18)", 
        "(0.18 < x < 0.21)"
    };

    // Title to add to each subplot
    std::string plotTitle = "0.16<z<0.22, 0.325<P_{T}<0.475, -0.5<x_{F}<-0.4, M_{x}>0.95";

    // Updated loop to accommodate the new plot order
    for (size_t i = 0; i < suffixes.size(); ++i) {
        c->cd(i + 1);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.15);

        bool firstGraphDrawn = false;

        // Create a legend in the top right with a border and background
        TLegend *legend = new TLegend(0.55, 0.7, 0.9, 0.9);  // Adjust position and size
        legend->SetTextSize(0.035);  // Adjust text size
        legend->SetBorderSize(1);  // Set border size
        legend->SetFillStyle(1001);   // Solid white background

        for (size_t p = 0; p < prefixes.size(); ++p) {
            std::string key = prefixes[p] + "chi2Fits" + suffixes[i];
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

                TGraphErrors *graphStat = createTGraphErrors(x, y, yStatErr, markers[p], 1.0, colors[p]);  // Adjusted marker size to 1.0
                graphStat->SetLineWidth(1);  // Ensure the line is centered properly
                
                setAxisLabelsAndRanges(graphStat, xLabel, yLabels[i], xLimits, 
                                       (suffixes[i] == "ALL") ? std::make_pair(-0.1, 0.5) :
                                       (suffixes[i] == "doubleratio") ? std::make_pair(-0.02, 0.4) :
                                       std::make_pair(-0.06, 0.06));  // Adjusted y-axis range for doubleratio
                
                if (!firstGraphDrawn) {
                    graphStat->Draw("AP");
                    firstGraphDrawn = true;
                } else {
                    graphStat->Draw("P SAME");
                }

                // Add each entry to the legend with the corresponding color and label
                TLegendEntry* legendEntry = legend->AddEntry(graphStat, legendLabels[p].c_str(), "p");
                legendEntry->SetTextColor(colors[p]);  // Set the text color to match the graph color
            }
        }

        // Add the dashed gray line at y = 0
        TLine *line = new TLine(xLimits.first, 0, xLimits.second, 0);
        line->SetLineColor(kGray+2);
        line->SetLineStyle(7);  // Dashed line
        line->Draw("same");

        // Draw the legend in the current subplot
        legend->Draw();

        // Add a title to each subplot
        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.035);
        latex->SetTextAlign(22);  // Centered text alignment
        latex->DrawLatex(0.535, 0.93, plotTitle.c_str());  // Adjust position and size
    }

    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs(outputFileName.c_str());
    delete c;
}

void plotNormalizedFLLOverFUU(
    const std::map<std::string, std::vector<std::vector<double>>> &asymmetryData,
    const std::map<std::string, std::vector<std::vector<double>>> &kinematicData,
    const std::string &outputFileName) {
    // Disable the display of titles globally
    gStyle->SetOptTitle(0);

    // Step 1: Fit the x-dependence of FLL/FUU
    auto it = asymmetryData.find("xchi2FitsALL");
    if (it == asymmetryData.end()) {
        std::cerr << "Error: xchi2FitsALL not found in asymmetry data.\n";
        return;
    }
    const auto &xData = it->second;

    std::vector<double> xVals, yVals, yErrs;
    for (const auto &entry : xData) {
        xVals.push_back(entry[0]); // x value
        yVals.push_back(entry[1]); // FLL/FUU value
        yErrs.push_back(entry[2]); // Error in FLL/FUU
    }

    TGraphErrors *graph = new TGraphErrors(xVals.size(), &xVals[0], &yVals[0], nullptr, &yErrs[0]);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.8);
    graph->SetMarkerColor(kBlack);
    graph->SetTitle(""); // Remove the title

    // Define the function y = x^a
    TF1 *fitFunc = new TF1("fitFunc", "pow(x, [0])", 0.06, 0.6);
    fitFunc->SetParameter(0, 1); // Initial guess for the exponent [0]
    fitFunc->SetParName(0, "Exponent");

    // Perform the fit
    graph->Fit(fitFunc, "Q"); // "Q" for quiet mode

    // Get fit parameter
    double exponent = fitFunc->GetParameter(0);

    // Step 2: Normalize FLL/FUU for other variables
    std::vector<std::string> variables = {"Mx", "z", "PT", "xF", "t"};
    std::map<std::string, TGraphErrors*> normalizedGraphs;

    for (const auto &var : variables) {
        std::string key = var + "chi2FitsALL";
        auto asymIt = asymmetryData.find(key);
        if (asymIt == asymmetryData.end()) {
            std::cerr << "Warning: " << key << " not found in asymmetry data.\n";
            continue;
        }
        const auto &asymData = asymIt->second;

        // Get corresponding kinematic data
        std::string kinKey = var + "Kinematics";
        auto kinIt = kinematicData.find(kinKey);
        if (kinIt == kinematicData.end()) {
            std::cerr << "Warning: " << kinKey << " not found in kinematic data.\n";
            continue;
        }
        const auto &kinData = kinIt->second;

        if (asymData.size() != kinData.size()) {
            std::cerr << "Warning: Size mismatch between asymmetry and kinematic data for " << var << ".\n";
            continue;
        }

        std::vector<double> varVals, normYVals, normYErrs;
        for (size_t i = 0; i < asymData.size(); ++i) {
            double varValue = asymData[i][0]; // Mean value of the variable
            double xValue = kinData[i][2];    // x value from kinematic data
            double yValue = asymData[i][1];   // FLL/FUU value
            double yErr = asymData[i][2];     // Error in FLL/FUU

            // For "t", we need to plot -t
            if (var == "t") {
                varValue = -varValue; // Convert t to -t
            }

            // Evaluate the fitted function at xValue
            double fittedY = fitFunc->Eval(xValue);

            // Normalize yValue
            double normY = yValue / fittedY;
            double normYErr = yErr / fabs(fittedY); // Error propagation assuming fittedY has negligible error

            varVals.push_back(varValue);
            normYVals.push_back(normY);
            normYErrs.push_back(normYErr);
        }

        // Create TGraphErrors for normalized data
        TGraphErrors *normGraph = new TGraphErrors(varVals.size(), &varVals[0], &normYVals[0], nullptr, &normYErrs[0]);
        normGraph->SetMarkerStyle(20);
        normGraph->SetMarkerSize(0.8);
        normGraph->SetMarkerColor(kBlack);
        normGraph->SetTitle(""); // Remove the title

        normalizedGraphs[var] = normGraph;
    }

    // Step 3: Plot the results
    TCanvas *c = new TCanvas("c", "Normalized FLL/FUU Plots", 1200, 800);
    c->Divide(3, 2);

    // Plot 1: FLL/FUU vs x with fit
    c->cd(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);

    setAxisLabelsAndRanges(graph, "x_{B}", "F_{LL}/F_{UU}", {0.06, 0.6}, {0, 0.6});
    graph->Draw("AP");

    // Draw the fit function
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("same");

    // Add legend in the top-left corner
    TLegend *leg = new TLegend(0.2, 0.75, 0.675, 0.9);
    leg->AddEntry(graph, "Data", "P");
    leg->AddEntry(fitFunc, Form("Fit: y = x^{%.3f}", exponent), "L");
    leg->Draw();

    // Draw the dashed gray line at y = 1
    double xmin_main = graph->GetXaxis()->GetXmin();
    double xmax_main = graph->GetXaxis()->GetXmax();
    TLine *line_main = new TLine(xmin_main, 1, xmax_main, 1);
    line_main->SetLineColor(kGray + 2);
    line_main->SetLineStyle(7); // Dashed line
    line_main->Draw("same");

    // Add text "M_{x} > 1.35 GeV" in the x plot
    TLatex *latex_main = new TLatex();
    latex_main->SetNDC();
    latex_main->SetTextSize(0.04);
    latex_main->SetTextAlign(33); // Bottom right alignment
    latex_main->DrawLatex(0.875, 0.225, "M_{x} > 1.35 GeV"); // Moved up slightly

    // Plot 2-6: Normalized FLL/FUU vs other variables
    int pad = 2;
    for (const auto &var : variables) {
        if (normalizedGraphs.find(var) == normalizedGraphs.end()) {
            continue;
        }

        c->cd(pad);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.15);

        TGraphErrors *normGraph = normalizedGraphs[var];

        // Set axis labels and ranges based on variable
        if (var == "Mx") {
            setAxisLabelsAndRanges(normGraph, "M_{x} (GeV)", "(F_{LL}/F_{UU}) / x^{a}", {0.0, 2.5}, {0.5, 1.2});
        } else if (var == "z") {
            setAxisLabelsAndRanges(normGraph, "z", "(F_{LL}/F_{UU}) / x^{a}", {0.0, 0.8}, {0.5, 1.2});
        } else if (var == "PT") {
            setAxisLabelsAndRanges(normGraph, "P_{T} (GeV)", "(F_{LL}/F_{UU}) / x^{a}", {0.0, 1.0}, {0.5, 1.2});
        } else if (var == "xF") {
            setAxisLabelsAndRanges(normGraph, "x_{F}", "(F_{LL}/F_{UU}) / x^{a}", {-0.8, 0.6}, {0.5, 1.2});
        } else if (var == "t") {
            setAxisLabelsAndRanges(normGraph, "-t (GeV^{2})", "(F_{LL}/F_{UU}) / x^{a}", {0.0, 8.0}, {0.5, 1.2});
        }

        normGraph->Draw("AP");

        // Draw the dashed gray line at y = 1
        double xmin = normGraph->GetXaxis()->GetXmin();
        double xmax = normGraph->GetXaxis()->GetXmax();
        TLine *line = new TLine(xmin, 1, xmax, 1);
        line->SetLineColor(kGray + 2);
        line->SetLineStyle(7); // Dashed line
        line->Draw("same");

        // Add text "M_{x} > 1.35 GeV" in the x, z, PT, xF, t plots but NOT in Mx plot
        if (var != "Mx") {
            TLatex *latex = new TLatex();
            latex->SetNDC();
            latex->SetTextSize(0.04);
            latex->SetTextAlign(33); // Bottom right alignment
            latex->DrawLatex(0.875, 0.225, "M_{x} > 1.35 GeV"); // Moved up slightly
        }

        pad++;
    }

    // Save the canvas
    gSystem->Exec("mkdir -p output/epX_plots");
    c->SaveAs(outputFileName.c_str());

    // Clean up
    delete c;
    delete fitFunc;
    delete graph;
    for (auto &entry : normalizedGraphs) {
        delete entry.second;
    }
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
    // plotDependence(asymmetryData, "Mx", "M_{x} (GeV)", {0.0, 2.5}, "output/epX_plots/Mx_dependence_plots.png");
    // plotDependence(asymmetryData, "x", "x_{B}", {0.06, 0.6}, "output/epX_plots/x_dependence_plots.png");
    // plotDependence(asymmetryData, "x", "x_{B}", {0.06, 0.6}, "output/epX_plots/x_dependence_plots_comparison.png", "xall");
    // plotDependence(asymmetryData, "z", "z", {0.0, 0.8}, "output/epX_plots/z_dependence_plots.png");
    // plotDependence(asymmetryData, "z", "z", {0.0, 0.8}, "output/epX_plots/z_dependence_plots_comparison.png", "zall");
    // plotDependence(asymmetryData, "PT", "P_{T} (GeV)", {0.0, 1.0}, "output/epX_plots/PT_dependence_plots.png");
    // plotDependence(asymmetryData, "PT", "P_{T} (GeV)", {0.0, 1.0}, "output/epX_plots/PT_dependence_plots_comparison.png", "PTall");
    // plotDependence(asymmetryData, "xF", "x_{F}", {-0.8, 0.6}, "output/epX_plots/xF_dependence_plots.png");
    // plotDependence(asymmetryData, "xF", "x_{F}", {-0.8, 0.6}, "output/epX_plots/xF_dependence_plots_comparison.png", "xFall");
    // plotDependence(asymmetryData, "Q2multi1", "Q^{2} (GeV^{2})", {1, 3.5}, "output/epX_plots/Q2multi1_dependence_plots.png");
    // plotDependence(asymmetryData, "Q2multi2", "Q^{2} (GeV^{2})", {1, 3.5}, "output/epX_plots/Q2multi2_dependence_plots.png");
    // plotDependence(asymmetryData, "Q2multi3", "Q^{2} (GeV^{2})", {1, 3.5}, "output/epX_plots/Q2multi3_dependence_plots.png");
    // plotDependence(asymmetryData, "Q2extramulti1", "Q^{2} (GeV^{2})", {1, 3.5}, "output/epX_plots/Q2extramulti1_dependence_plots.png");
    // plotDependence(asymmetryData, "Q2extramulti2", "Q^{2} (GeV^{2})", {1, 3.5}, "output/epX_plots/Q2extramulti2_dependence_plots.png");
    // plotDependence(asymmetryData, "Q2extramulti3", "Q^{2} (GeV^{2})", {1, 3.5}, "output/epX_plots/Q2extramulti3_dependence_plots.png");
    plotDependence(asymmetryData, "t", "t (GeV^{2})", {-8, 0}, "output/epX_plots/t_dependence_plots.png");

    // Call the new function
    plotNormalizedFLLOverFUU(asymmetryData, kinematicData, "output/epX_plots/normalized_FLL_over_FUU.png");

    // std::vector<std::string> prefixes = {"Q2multi1", "Q2multi2", "Q2multi3"};
    // plotMultipleQ2multiDependence(asymmetryData, prefixes, "Q^{2} (GeV^{2})", {1, 3.5}, "output/epX_plots/Q2multi_dependence_plots.png");

    // std::vector<std::string> extraprefixes = {"Q2extramulti1", "Q2extramulti2", "Q2extramulti3"};
    // plotMultipleQ2extramultiDependence(asymmetryData, extraprefixes, "Q^{2} (GeV^{2})", {1, 3.5}, "output/epX_plots/Q2extramulti_dependence_plots.png");

    // plotRunnumDependence(asymmetryData, "runnum", "run number", "output/epX_plots/runnum_dependence_plots.png");

    // plotTargetPolarizationDependence(targetPolarizationData, "run number", "output/epX_plots/target_polarization_dependence.png");

    return 0;
}
