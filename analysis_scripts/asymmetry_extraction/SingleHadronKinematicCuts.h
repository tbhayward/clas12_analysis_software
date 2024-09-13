#pragma once

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTree.h>
#include <optional>

#include "BaseKinematicCuts.h"

class SingleHadronKinematicCuts : public BaseKinematicCuts {
public:
    // Constructor that takes TTreeReader and TTree to check for branch existence
    SingleHadronKinematicCuts(TTreeReader& reader, TTree* tree);

    // Overridden function for applying cuts
    bool applyCuts(int currentFits, bool isMC) override;

    // Destructor
    ~SingleHadronKinematicCuts() = default;

private:
    // Required variables
    TTreeReaderValue<int> runnum;
    TTreeReaderValue<double> e_theta;
    TTreeReaderValue<double> e_phi;
    TTreeReaderValue<double> vz_e;
    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> Mx;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> z;
    TTreeReaderValue<double> pT;
    TTreeReaderValue<double> xF;
    TTreeReaderValue<double> target_pol;

    // Optional variables (non-null only if the corresponding branches exist)
    std::optional<TTreeReaderValue<double>> Mx1;
    std::optional<TTreeReaderValue<double>> Mx2;
    std::optional<TTreeReaderValue<double>> Mx23;
};