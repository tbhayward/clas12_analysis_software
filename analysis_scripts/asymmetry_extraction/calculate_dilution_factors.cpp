#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <TLatex.h>
#include <string>
#include <iostream>

// tbhayward libraries
#include "common_vars.h"
#include "load_bins_from_csv.h"
#include "load_run_info_from_csv.h"
#include "dilution_factor.h"
#include "calculate_dilution_factors.h"
#include "asymmetry_fits.h"
#include "BaseKinematicCuts.h"
#include "KinematicCuts.h"
#include "InclusiveKinematicCuts.h"
#include "SingleHadronKinematicCuts.h"
#include "B2BDihadronKinematicCuts.h"
#include "DihadronKinematicCuts.h"
#include "dvcsKinematicCuts.h"
#include "eppi0KinematicCuts.h"
#include "GeneralExclusiveKinematicCuts.h"
#include "formatLabelName.h"
#include "readChi2Fits.h"
#include "histConfigs.h"
#include "charge_accumulation.h"
#include "plot_data.h"
#include "modifyTree.h"
#include "fitting_process.h"

// Select dataset: 1 = RGC Su22, 2 = RGC Fa22, 3 = RGC Sp23
constexpr int data_set = 2;

struct DataSetConfig {
    std::string name;
    double xA, xC, xCH, xHe, xf;
    std::string nh3File;
    std::string cFile;
    std::string chFile;
    std::string heFile;
    std::string emptyFile;
};

const std::vector<DataSetConfig> dataSetConfigs = {
    {
        "RGC_Su22",
        0.7624, 0.0747, 0.0390, 0.0249, 0.0989,
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_NH3_eX.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_C_eX.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_CH2_eX.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_He_eX.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_ET_eX.root"
    },
    {
        "RGC_Fa22",
        0.5839, 0.1992, 0.1802, 0.0301, 0.0066,
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_C_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_CH2_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_He_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_ET_epi+_2.root"
    },
    {
        "RGC_Sp23",
        0.4543, 0.1509, 0.1966, 0.1205, 0.0777,
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_NH3_eX.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_C_eX.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_CH2_eX.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_He_eX.root",
        "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_ET_eX.root"
    },
    {
        "RGC_Su22_enpi+",
        0.7624, 0.0747, 0.0390, 0.0249, 0.0989,
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_NH3_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_C_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_CH2_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_He_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_ET_epi+_2.root"
    },
    {
        "RGC_Fa22_enpi+",
        0.5839, 0.1992, 0.1802, 0.0301, 0.0066,
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_C_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_CH2_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_He_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_ET_epi+_2.root"
    },
    {
        "RGC_Sp23_enpi+",
        0.4543, 0.1509, 0.1966, 0.1205, 0.0777,
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_NH3_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_C_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_CH2_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_He_epi+_2.root",
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_ET_epi+_2.root"
    },
    {
        "RGC_Fa22_enpi+_ISR_FSR_0",
        0.5839, 0.1992, 0.1802, 0.0301, 0.0066,
        "/work/clas12/thayward/ISR_FSR_study/rgc_fa22_inb_NH3_epi+_ISR_FSR_0_2.root",
        "/work/clas12/thayward/ISR_FSR_study/rgc_fa22_inb_C_epi+_ISR_FSR_0_2.root",
        "/work/clas12/thayward/ISR_FSR_study/rgc_fa22_inb_CH2_epi+_ISR_FSR_0_2.root",
        "/work/clas12/thayward/ISR_FSR_study/rgc_fa22_inb_He_epi+_ISR_FSR_0_2.root",
        "/work/clas12/thayward/ISR_FSR_study/rgc_fa22_inb_ET_epi+_ISR_FSR_0_2.root",
    },
};

// NH3 periods defined as a pair of run numbers (start, end)
static const std::vector<std::pair<int, int>> nh3_periods = {
    {16137, 16148}, {16156, 16178}, {16211, 16228}, {16231, 16260},
    {16318, 16333}, {16335, 16357}, {16709, 16720}, {16721, 16766}, {16767, 16772}
};

// Function to calculate standard deviation of a vector
double calculate_standard_deviation(const std::vector<double>& values) {
    if (values.size() < 2) return 0.0;
    double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double sq_sum = std::inner_product(
        values.begin(), values.end(),
        values.begin(), 0.0,
        std::plus<double>(),
        [&](double /*unused*/, double b) { return (b - mean) * (b - mean); }
    );
    return std::sqrt(sq_sum / (values.size() - 1));
}

double calculate_dilution_error(double nA, double nC, double nCH, double nMT, double nf,
                                double xA, double xC, double xCH, double xHe, double xf) {
    // // First part of the expression (term1)
    // double term1 = 5438.81 * nA * nf * std::pow(xC, 2) * std::pow(xCH, 2)
    //                * std::pow(xf, 2) * std::pow(xHe, 2)
    //                * std::pow((-1.0 * nMT * xA + nA * xHe), 2)
    //                * std::pow((1.0 * nMT * xC * xCH
    //                            - 1.19655 * nCH * xC * xHe
    //                            + 0.196547 * nC * xCH * xHe), 2);

    // // Second part of the expression (term2)
    // double term2 = 85044.9 * nA * nCH * std::pow(xC, 2) * std::pow(xCH, 2)
    //                * std::pow(xf, 2) * std::pow(xHe, 2)
    //                * std::pow((-1.0 * nMT * xA + nA * xHe), 2)
    //                * std::pow((1.0 * nMT * xC * xf
    //                            - 0.302592 * nf * xC * xHe
    //                            - 0.697408 * nC * xf * xHe), 2);

    // // Third part of the expression (term3)
    // double term3 = 47470.0 * nA * nC * std::pow(xC, 2) * std::pow(xCH, 2)
    //                * std::pow(xf, 2) * std::pow(xHe, 2)
    //                * std::pow((-1.0 * nMT * xA + nA * xHe), 2)
    //                * std::pow((1.0 * nMT * xCH * xf
    //                            - 0.0665285 * nf * xCH * xHe
    //                            - 0.933472 * nCH * xf * xHe), 2);

    // // Fourth part of the expression (term4)
    // double term4 = 1371.83 * std::pow(nMT, 2) * std::pow(xA, 2)
    //                * std::pow((1.0 * nMT * xC * xCH * xf
    //                            + (-1.97748 * nf * xC * xCH
    //                               + 6.62306 * nCH * xC * xf
    //                               - 5.64558 * nC * xCH * xf) * xHe), 2)
    //                * std::pow((1.0 * nMT * xC * xCH * xf
    //                            + (0.0136533 * nf * xC * xCH
    //                               - 1.25054 * nCH * xC * xf
    //                               + 0.236883 * nC * xCH * xf) * xHe), 2);

    // // Fifth part of the expression (term5)
    // double term5_numerator_part1 = 73.2426 * std::pow(nMT, 2) * xA
    //                                * std::pow(xC, 2) * std::pow(xCH, 2)
    //                                * std::pow(xf, 2);

    // double term5_numerator_part2 = nMT * xC * xCH * xf
    //                                * (2.0 * nf * xA * xC * xCH
    //                                   + (-183.185 * nCH * xA * xC
    //                                      + 34.6998 * nC * xA * xCH
    //                                      + 1.42109e-14 * nA * xC * xCH) * xf);

    // double term5_numerator_part3 = (-1.97748 * std::pow(nf, 2) * xA
    //                                * std::pow(xC, 2) * std::pow(xCH, 2))
    //                                + nf * xC * xCH
    //                                  * (187.746 * nCH * xA * xC
    //                                     - 39.9547 * nC * xA * xCH
    //                                     - 145.836 * nA * xC * xCH) * xf
    //                                + (-606.623 * std::pow(nCH, 2) * xA
    //                                   * std::pow(xC, 2)
    //                                   + nCH * xC
    //                                     * (632.002 * nC * xA
    //                                        + 576.683 * nA * xC) * xCH
    //                                   + nC
    //                                     * (-97.9502 * nC * xA
    //                                        - 430.847 * nA * xC)
    //                                     * std::pow(xCH, 2)) * std::pow(xf, 2);

    // double term5_numerator = term5_numerator_part1
    //                          + term5_numerator_part2 * xHe
    //                          + term5_numerator_part3 * std::pow(xHe, 2);

    // double term5 = 0.255725 * nA * nMT * std::pow(term5_numerator, 2);

    // double denominator = std::pow(nA, 3) * std::pow(xHe, 2)
    //                      * std::pow((73.2426 * nMT * xC * xCH * xf
    //                                  + 1.0 * nf * xC * xCH * xHe
    //                                  - 91.5925 * nCH * xC * xf * xHe
    //                                  + 17.3499 * nC * xCH * xf * xHe), 4);

    // double sigma_df = 27.3473
    //                   * std::sqrt((term1 + term2 + term3 + term4 + term5) / denominator);
    // return sigma_df;
    // pre thermal contraction calculations

    // Precompute common factors
    double commonFactor = -1.0 * nMT * xA + nA * xHe;
    double xC2 = std::pow(xC, 2);
    double xCH2 = std::pow(xCH, 2);
    double xf2 = std::pow(xf, 2);
    double xHe2 = std::pow(xHe, 2);

    // Term 1: Remove erroneous 12.3729 factor
    double factor1b = nMT * xC * xCH - 1.32754 * nCH * xC * xHe + 0.327539 * nC * xCH * xHe;
    double term1 = 1246.13 * nA * nf * xC2 * xCH2 * xf2 * xHe2 * 
                   std::pow(commonFactor, 2) * std::pow(factor1b, 2);

    // Term 2 (unchanged)
    double factor2b = nMT * xC * xf - 0.31065 * nf * xC * xHe - 0.68935 * nC * xf * xHe;
    double term2 = 22757.0 * nA * nCH * xC2 * xCH2 * xf2 * xHe2 * 
                   std::pow(commonFactor, 2) * std::pow(factor2b, 2);

    // Term 3 (unchanged)
    double factor3b = nMT * xCH * xf - 0.10006 * nf * xCH * xHe - 0.89994 * nCH * xf * xHe;
    double term3 = 13352.7 * nA * nC * xC2 * xCH2 * xf2 * xHe2 * 
                   std::pow(commonFactor, 2) * std::pow(factor3b, 2);

    // Term 4 (unchanged)
    double factor4a = nMT * xC * xCH * xf + 
                     xHe * (-1.72606 * nf * xC * xCH + 6.04866 * nCH * xC * xf - 5.3226 * nC * xCH * xf);
    double factor4b = nMT * xC * xCH * xf + 
                     xHe * (-0.0278707 * nf * xC * xCH - 1.20844 * nCH * xC * xf + 0.236306 * nC * xCH * xf);
    double term4 = 432.106 * std::pow(nMT, 2) * std::pow(xA, 2) * 
                   std::pow(factor4a, 2) * std::pow(factor4b, 2);

    // Term 5 (unchanged)
    double inner5 = 
        35.88 * std::pow(nMT, 2) * xA * xC2 * xCH2 * xf2 +
        nMT * xA * xC * xCH * xf * xHe * 
            (-2.0 * nf * xC * xCH - 86.7173 * nCH * xC * xf + 16.9573 * nC * xCH * xf) +
        xHe2 * (
            1.72606 * std::pow(nf, 2) * xA * xC2 * xCH2 +
            nf * xC * xCH * xf * 
                (68.7911 * nCH * xA * xC - 9.31212 * nC * xA * xCH - 60.9311 * nA * xC * xCH) +
            xf2 * (
                -262.262 * std::pow(nCH, 2) * xA * xC2 +
                nCH * xC * xCH * (282.065 * nC * xA + 260.384 * nA * xC) +
                nC * xCH2 * (-45.1285 * nC * xA - 199.453 * nA * xC)
            )
        );
    double term5 = 0.33565 * nA * nMT * std::pow(inner5, 2);

    // Combine numerator terms
    double numerator = term1 + term2 + term3 + term4 + term5;

    // Compute denominator polynomial
    double poly = 
        35.88 * nMT * xC * xCH * xf -
        1.0 * nf * xC * xCH * xHe -
        43.3586 * nCH * xC * xf * xHe +
        8.47866 * nC * xCH * xf * xHe;

    // Handle potential division by zero
    if (std::abs(poly) < 1e-10) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Compute full denominator
    double denominator = std::pow(nA, 3) * xHe2 * std::pow(poly, 4);

    // Apply global multiplicative factor 12.3729 AFTER sqrt
    double result = std::sqrt(numerator / denominator);
    return 12.3729 * result;
}


std::pair<double, double> calculate_dilution_and_error(
    double nA, double nC, double nCH, double nMT, double nf,
    double xA, double xC, double xCH, double xHe, double xf) {
    // double dilution = (27.3473 * (-1.0 * nMT * xA + nA * xHe)
    //                    * (-0.505693 * nMT * xC * xCH * xf
    //                       + (1.0 * nf * xC * xCH
    //                          - 3.34924 * nCH * xC * xf
    //                          + 2.85493 * nC * xCH * xf) * xHe))
    //                   / (nA * xHe
    //                      * (73.2426 * nMT * xC * xCH * xf
    //                         + 1.0 * nf * xC * xCH * xHe
    //                         - 91.5925 * nCH * xC * xf * xHe
    //                         + 17.3499 * nC * xCH * xf * xHe));
    // the above is the calculation for the old numbers before the thermal contraction was 
    // considered

    double dilution = (12.3729 * (-1.0 * nMT * xA + nA * xHe)
                       * (-0.579353 * nMT * xC * xCH * xf
                          + (1.0 * nf * xC * xCH
                             - 3.50431 * nCH * xC * xf
                             + 3.08366 * nC * xCH * xf) * xHe))
                      / (nA * xHe
                         * (35.88 * nMT * xC * xCH * xf
                            - 1.0 * nf * xC * xCH * xHe
                            - 43.3586 * nCH * xC * xf * xHe
                            + 8.47866 * nC * xCH * xf * xHe));

    // double packing_fraction = (0.699832)*(nA/xA - nMT/xHe)
    //     / (1.25055*nCH/xCH - 0.23688*nC/xC - 0.013668*nf/xf - nMT/xHe);

    std::cout << dilution << " ";

    double error = calculate_dilution_error(nA, nC, nCH, nMT, nf, xA, xC, xCH, xHe, xf);
    return { dilution, error };
}

std::vector<std::pair<double, double>> calculate_dilution_factors() {
    if (data_set < 1 || data_set > static_cast<int>(dataSetConfigs.size())) {
        std::cerr << "ERROR: data_set=" << data_set
                  << " is invalid; must be 1–" << dataSetConfigs.size() << "\n";
        return {};
    }
    const auto& dsCfg = dataSetConfigs[data_set - 1];

    // Load ROOT files and trees
    TFile* nh3File   = TFile::Open(dsCfg.nh3File.c_str());
    TFile* cFile     = TFile::Open(dsCfg.cFile.c_str());
    TFile* chFile    = TFile::Open(dsCfg.chFile.c_str());
    TFile* heFile    = TFile::Open(dsCfg.heFile.c_str());
    TFile* emptyFile = TFile::Open(dsCfg.emptyFile.c_str());

    TTree* nh3       = static_cast<TTree*>(nh3File->Get("PhysicsEvents"));
    TTree* treeC     = static_cast<TTree*>(cFile->Get("PhysicsEvents"));
    TTree* treeCh    = static_cast<TTree*>(chFile->Get("PhysicsEvents"));
    TTree* treeHe    = static_cast<TTree*>(heFile->Get("PhysicsEvents"));
    TTree* treeEmpty = static_cast<TTree*>(emptyFile->Get("PhysicsEvents"));

    TTreeReader nh3Reader(nh3), cReader(treeC), chReader(treeCh),
                  heReader(treeHe), emptyReader(treeEmpty);
    TTreeReaderValue<int> runnum(nh3Reader, "runnum");

    // Allocate kinematic cuts
    BaseKinematicCuts* nh3Cuts   = nullptr;
    BaseKinematicCuts* cCuts     = nullptr;
    BaseKinematicCuts* chCuts    = nullptr;
    BaseKinematicCuts* heCuts    = nullptr;
    BaseKinematicCuts* emptyCuts = nullptr;
    switch (channel) {
        case 0:
            nh3Cuts   = new InclusiveKinematicCuts(nh3Reader);
            cCuts     = new InclusiveKinematicCuts(cReader);
            chCuts    = new InclusiveKinematicCuts(chReader);
            heCuts    = new InclusiveKinematicCuts(heReader);
            emptyCuts = new InclusiveKinematicCuts(emptyReader);
            break;
        case 1:
            nh3Cuts   = new SingleHadronKinematicCuts(nh3Reader);
            cCuts     = new SingleHadronKinematicCuts(cReader);
            chCuts    = new SingleHadronKinematicCuts(chReader);
            heCuts    = new SingleHadronKinematicCuts(heReader);
            emptyCuts = new SingleHadronKinematicCuts(emptyReader);
            break;
        case 2:
            nh3Cuts   = new B2BDihadronKinematicCuts(nh3Reader);
            cCuts     = new B2BDihadronKinematicCuts(cReader);
            chCuts    = new B2BDihadronKinematicCuts(chReader);
            heCuts    = new B2BDihadronKinematicCuts(heReader);
            emptyCuts = new B2BDihadronKinematicCuts(emptyReader);
            break;
        case 3:
            nh3Cuts   = new DihadronKinematicCuts(nh3Reader);
            cCuts     = new DihadronKinematicCuts(cReader);
            chCuts    = new DihadronKinematicCuts(chReader);
            heCuts    = new DihadronKinematicCuts(heReader);
            emptyCuts = new DihadronKinematicCuts(emptyReader);
            break;
        case 4:
            nh3Cuts   = new dvcsKinematicCuts(nh3Reader);
            cCuts     = new dvcsKinematicCuts(cReader);
            chCuts    = new dvcsKinematicCuts(chReader);
            heCuts    = new dvcsKinematicCuts(heReader);
            emptyCuts = new dvcsKinematicCuts(emptyReader);
            break;
        case 5:
            nh3Cuts   = new eppi0KinematicCuts(nh3Reader);
            cCuts     = new eppi0KinematicCuts(cReader);
            chCuts    = new eppi0KinematicCuts(chReader);
            heCuts    = new eppi0KinematicCuts(heReader);
            emptyCuts = new eppi0KinematicCuts(emptyReader);
            break;
        case 6:
            nh3Cuts   = new GeneralExclusiveKinematicCuts(nh3Reader);
            cCuts     = new GeneralExclusiveKinematicCuts(cReader);
            chCuts    = new GeneralExclusiveKinematicCuts(chReader);
            heCuts    = new GeneralExclusiveKinematicCuts(heReader);
            emptyCuts = new GeneralExclusiveKinematicCuts(emptyReader);
            break;
        default:
            std::cerr << "Invalid channel specified." << std::endl;
            return {};
    }

    std::vector<std::pair<double, double>> dilutionResults;
    TGraphErrors* gr_dilution = new TGraphErrors();

    auto fill_hist = [&](TTreeReader& reader, TH1D* hist, BaseKinematicCuts* cuts) {
        TTreeReaderValue<double> var(reader, propertyNames[currentFits].c_str());
        while (reader.Next()) {
            if (!cuts->applyCuts(currentFits, false)) continue;
            hist->Fill(*var);
        }
        reader.Restart();
    };

    size_t nBins = allBins[currentFits].size() - 1;
    for (size_t i = 0; i < nBins; ++i) {
        double varMin = allBins[currentFits][i];
        double varMax = allBins[currentFits][i + 1];

        TH1D* h_nh3   = new TH1D("h_nh3",    "", 1, varMin, varMax);
        TH1D* h_c     = new TH1D("h_c",      "", 1, varMin, varMax);
        TH1D* h_ch    = new TH1D("h_ch",     "", 1, varMin, varMax);
        TH1D* h_he    = new TH1D("h_he",     "", 1, varMin, varMax);
        TH1D* h_empty = new TH1D("h_empty",  "", 1, varMin, varMax);

        fill_hist(nh3Reader,   h_nh3,   nh3Cuts);
        fill_hist(cReader,     h_c,     cCuts);
        fill_hist(chReader,    h_ch,    chCuts);
        fill_hist(heReader,    h_he,    heCuts);
        fill_hist(emptyReader, h_empty, emptyCuts);

        double nA  = h_nh3->GetBinContent(1);
        double nC  = h_c->GetBinContent(1);
        double nCH = h_ch->GetBinContent(1);
        double nMT = h_he->GetBinContent(1);
        double nf  = h_empty->GetBinContent(1);

        // only for Sp23, empty target was 20% He → scale its yield
        if (data_set == 3) {
            constexpr double fa22_ET = 0.602465;
            constexpr double sp23_ET = 0.733067;
            nf *= (fa22_ET / sp23_ET);
        }

        auto [dil, err] = calculate_dilution_and_error(
            nA, nC, nCH, nMT, nf,
            dsCfg.xA, dsCfg.xC, dsCfg.xCH, dsCfg.xHe, dsCfg.xf
        );

        dilutionResults.emplace_back(dil, err);
        gr_dilution->SetPoint(i, (varMin + varMax) / 2.0, dil);
        gr_dilution->SetPointError(i, 0, err);

        delete h_nh3; delete h_c; delete h_ch; delete h_he; delete h_empty;
    }

    TCanvas* canvasDF = new TCanvas("c_dilution", "Dilution Factor Plot", 800, 600);
    canvasDF->SetLeftMargin(0.15);
    canvasDF->SetBottomMargin(0.15);

    std::string prefix = propertyNames[currentFits];
    HistConfig cfg = histConfigs.count(prefix) ? histConfigs[prefix] : HistConfig{100, 0, 1};
    gr_dilution->SetTitle("");
    gr_dilution->GetXaxis()->SetTitle(formatLabelName(prefix).c_str());
    gr_dilution->GetXaxis()->SetLimits(cfg.xMin, cfg.xMax);
    gr_dilution->GetXaxis()->SetTitleSize(0.05);
    gr_dilution->GetYaxis()->SetTitle("D_{f}");
    gr_dilution->GetYaxis()->SetTitleSize(0.05);
    gr_dilution->GetYaxis()->SetTitleOffset(1.6);
    gr_dilution->GetYaxis()->SetRangeUser(0.10, 0.35);
    gr_dilution->SetMarkerStyle(20);
    gr_dilution->SetMarkerColor(kBlack);
    gr_dilution->Draw("AP");

    std::string outDir  = "output/dilution_factor_plots/";
    std::string outFile = outDir + "df_" + binNames[currentFits] + "_" + prefix + ".pdf";
    canvasDF->SaveAs(outFile.c_str());

    delete nh3Cuts; delete cCuts; delete chCuts; delete heCuts; delete emptyCuts;
    delete gr_dilution;
    delete canvasDF;

    return dilutionResults;
}