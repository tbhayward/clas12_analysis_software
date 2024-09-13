#pragma once

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTree.h>  // Include for the TTree usage

#include "BaseKinematicCuts.h"

class SingleHadronKinematicCuts : public BaseKinematicCuts {
public:
    // Constructor that takes TTreeReader and TTree to check for branch existence
    SingleHadronKinematicCuts(TTreeReader& reader, TTree* tree);

    // Overridden function for applying cuts
    bool applyCuts(int currentFits, bool isMC) override;

    // Destructor to clean up dynamically allocated memory
    ~SingleHadronKinematicCuts();

private:
    // Required variables
    TTreeReaderValue<int> runnum;
    TTreeReaderValue<double> e_theta;
    TTreeReaderValue<double> e_phi;
    TTreeReaderValue<double> vz_e;
    TTreeReaderValue<double> p_p;
    TTreeReaderValue<double> p_theta;
    TTreeReaderValue<double> p_phi;
    TTreeReaderValue<double> vz_p;
    TTreeReaderValue<double> Q2;
    TTreeReaderValue<double> W;
    TTreeReaderValue<double> Mx;
    TTreeReaderValue<double> x;
    TTreeReaderValue<double> y;
    TTreeReaderValue<double> z;
    TTreeReaderValue<double> pT;
    TTreeReaderValue<double> xF;
    TTreeReaderValue<double> phi;
    TTreeReaderValue<double> phi2;
    TTreeReaderValue<double> target_pol;

    // Optional variables (pointers initialized as nullptr)
    TTreeReaderValue<double>* Mx1 = nullptr;
    TTreeReaderValue<double>* Mx2 = nullptr;
    TTreeReaderValue<double>* Mx23 = nullptr;
};