#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <iostream>
#include <string>

int DetermineQ2yBin(float Q2, float y) {
    // Q2-y Bins 1-4
    if (Q2 > 1.000 && Q2 <= 2.000) {
        if (y > 0.650 && y <= 0.750) return 1;
        if (y > 0.550 && y <= 0.650) return 2;
        if (y > 0.450 && y <= 0.550) return 3;
        if (y > 0.300 && y <= 0.450) return 4;
    }
    // Q2-y Bins 5-8
    else if (Q2 > 2.000 && Q2 <= 3.000) {
        if (y > 0.650 && y <= 0.750) return 5;
        if (y > 0.550 && y <= 0.650) return 6;
        if (y > 0.450 && y <= 0.550) return 7;
        if (y > 0.300 && y <= 0.450) return 8;
    }
    // Q2-y Bins 9-12
    else if (Q2 > 3.000 && Q2 <= 4.000) {
        if (y > 0.650 && y <= 0.750) return 9;
        if (y > 0.550 && y <= 0.650) return 10;
        if (y > 0.450 && y <= 0.550) return 11;
        if (y > 0.300 && y <= 0.450) return 12;
    }
    // Q2-y Bins 13-15
    else if (Q2 > 4.000 && Q2 <= 5.000) {
        if (y > 0.650 && y <= 0.750) return 13;
        if (y > 0.550 && y <= 0.650) return 14;
        if (y > 0.450 && y <= 0.550) return 15;
    }
    // Q2-y Bins 16-17
    else if (Q2 > 5.000 && Q2 <= 7.000) {
        if (y > 0.650 && y <= 0.750) return 16;
        if (y > 0.550 && y <= 0.650) return 17;
    }
    return 0;
}

void process_file(const char* input_filename) {
    TFile *input_file = TFile::Open(input_filename);
    if (!input_file || input_file->IsZombie()) {
        std::cerr << "Error opening file: " << input_filename << std::endl;
        return;
    }

    TTree *input_tree;
    input_file->GetObject("PhysicsEvents", input_tree);
    if (!input_tree) {
        std::cerr << "Error: PhysicsEvents tree not found!" << std::endl;
        input_file->Close();
        return;
    }

    // Define the output filenames
    std::string base_name(input_filename);
    size_t dot_pos = base_name.find_last_of('.');
    if (dot_pos != std::string::npos) {
        base_name = base_name.substr(0, dot_pos);
    }

    // Declare variables for branches
    int runnum_val;
    int evnum_val;
    int helicity_val;
    double beam_pol_val;
    double target_pol_val;
    double e_p_val;
    double e_theta_val;
    double e_phi_val;
    double vz_e_val;
    double p_p_val;
    double p_theta_val;
    double p_phi_val;
    double vz_p_val;
    double Q2_val;
    double W_val;
    double Mx_val;
    double Mx2_val;
    double x_val;
    double y_val;
    double t_val;
    double tmin_val;
    double z_val;
    double xF_val;
    double pT_val;
    double zeta_val;
    double eta_val;
    double phi_val;
    double DepA_val;
    double DepB_val;
    double DepC_val;
    double DepV_val;
    double DepW_val;

    // Create output files and trees for each Q2-y bin
    std::vector<TFile*> output_files(18);
    std::vector<TTree*> output_trees(18);
    for (int i = 0; i < 18; ++i) {
        std::string output_filename = base_name + "_" + std::to_string(i) + ".root";
        output_files[i] = TFile::Open(output_filename.c_str(), "RECREATE");
        output_trees[i] = new TTree("PhysicsEvents", "PhysicsEvents");


        // Create branches for each variable
        output_trees[i]->Branch("runnum", &runnum_val, "runnum/I");
        output_trees[i]->Branch("evnum", &evnum_val, "evnum/I");
        output_trees[i]->Branch("helicity", &helicity_val, "helicity/I");
        output_trees[i]->Branch("beam_pol", &beam_pol_val, "beam_pol/D");
        output_trees[i]->Branch("target_pol", &target_pol_val, "target_pol/D");
        output_trees[i]->Branch("e_p", &e_p_val, "e_p/D");
        output_trees[i]->Branch("e_theta", &e_theta_val, "e_theta/D");
        output_trees[i]->Branch("e_phi", &e_phi_val, "e_phi/D");
        output_trees[i]->Branch("vz_e", &vz_e_val, "vz_e/D");
        output_trees[i]->Branch("p_p", &p_p_val, "p_p/D");
        output_trees[i]->Branch("p_theta", &p_theta_val, "p_theta/D");
        output_trees[i]->Branch("p_phi", &p_phi_val, "p_phi/D");
        output_trees[i]->Branch("vz_p", &vz_p_val, "vz_p/D");
        output_trees[i]->Branch("Q2", &Q2_val, "Q2/D");
        output_trees[i]->Branch("W", &W_val, "W/D");
        output_trees[i]->Branch("Mx", &Mx_val, "Mx/D");
        output_trees[i]->Branch("Mx2", &Mx2_val, "Mx2/D");
        output_trees[i]->Branch("x", &x_val, "x/D");
        output_trees[i]->Branch("y", &y_val, "y/D");
        output_trees[i]->Branch("t", &t_val, "t/D");
        output_trees[i]->Branch("tmin", &tmin_val, "tmin/D");
        output_trees[i]->Branch("z", &z_val, "z/D");
        output_trees[i]->Branch("xF", &xF_val, "xF/D");
        output_trees[i]->Branch("pT", &pT_val, "pT/D");
        output_trees[i]->Branch("zeta", &zeta_val, "zeta/D");
        output_trees[i]->Branch("eta", &eta_val, "eta/D");
        output_trees[i]->Branch("phi", &phi_val, "phi/D");
        output_trees[i]->Branch("DepA", &DepA_val, "DepA/D");
        output_trees[i]->Branch("DepB", &DepB_val, "DepB/D");
        output_trees[i]->Branch("DepC", &DepC_val, "DepC/D");
        output_trees[i]->Branch("DepV", &DepV_val, "DepV/D");
        output_trees[i]->Branch("DepW", &DepW_val, "DepW/D");
    }

    // Create a TTreeReader to read the input tree
    TTreeReader reader(input_tree);
    TTreeReaderValue<int> runnum(reader, "runnum");
    TTreeReaderValue<int> evnum(reader, "evnum");
    TTreeReaderValue<int> helicity(reader, "helicity");
    TTreeReaderValue<double> beam_pol(reader, "beam_pol");
    TTreeReaderValue<double> target_pol(reader, "target_pol");
    TTreeReaderValue<double> e_p(reader, "e_p");
    TTreeReaderValue<double> e_theta(reader, "e_theta");
    TTreeReaderValue<double> e_phi(reader, "e_phi");
    TTreeReaderValue<double> vz_e(reader, "vz_e");
    TTreeReaderValue<double> p_p(reader, "p_p");
    TTreeReaderValue<double> p_theta(reader, "p_theta");
    TTreeReaderValue<double> p_phi(reader, "p_phi");
    TTreeReaderValue<double> vz_p(reader, "vz_p");
    TTreeReaderValue<double> Q2(reader, "Q2");
    TTreeReaderValue<double> W(reader, "W");
    TTreeReaderValue<double> Mx(reader, "Mx");
    TTreeReaderValue<double> Mx2(reader, "Mx2");
    TTreeReaderValue<double> x(reader, "x");
    TTreeReaderValue<double> y(reader, "y");
    TTreeReaderValue<double> t(reader, "t");
    TTreeReaderValue<double> tmin(reader, "tmin");
    TTreeReaderValue<double> z(reader, "z");
    TTreeReaderValue<double> xF(reader, "xF");
    TTreeReaderValue<double> pT(reader, "pT");
    TTreeReaderValue<double> zeta(reader, "zeta");
    TTreeReaderValue<double> eta(reader, "eta");
    TTreeReaderValue<double> phi(reader, "phi");
    TTreeReaderValue<double> DepA(reader, "DepA");
    TTreeReaderValue<double> DepB(reader, "DepB");
    TTreeReaderValue<double> DepC(reader, "DepC");
    TTreeReaderValue<double> DepV(reader, "DepV");
    TTreeReaderValue<double> DepW(reader, "DepW");

    // Loop over entries and fill the corresponding output trees
    while (reader.Next()) {
        if (*Mx > 1.4) {
            // Assign values to variables
            int runnum_val = *runnum;
            int evnum_val = *evnum;
            int helicity_val = *helicity;
            double beam_pol_val = *beam_pol;
            double target_pol_val = *target_pol;
            double e_p_val = *e_p;
            double e_theta_val = *e_theta;
            double e_phi_val = *e_phi;
            double vz_e_val = *vz_e;
            double p_p_val = *p_p;
            double p_theta_val = *p_theta;
            double p_phi_val = *p_phi;
            double vz_p_val = *vz_p;
            double Q2_val = *Q2;
            double W_val = *W;
            double Mx_val = *Mx;
            double Mx2_val = *Mx2;
            double x_val = *x;
            double y_val = *y;
            double t_val = *t;
            double tmin_val = *tmin;
            double z_val = *z;
            double xF_val = *xF;
            double pT_val = *pT;
            double zeta_val = *zeta;
            double eta_val = *eta;
            double phi_val = *phi;
            double DepA_val = *DepA;
            double DepB_val = *DepB;
            double DepC_val = *DepC;
            double DepV_val = *DepV;
            double DepW_val = *DepW;

            // Determine the Q2-y bin and fill the corresponding tree
            int bin = DetermineQ2yBin(Q2_val, y_val);
            if (bin > 0 && bin < 18) {
                output_trees[bin]->Fill();
            }
            // Also fill the general tree for events passing Mx cut
            output_trees[0]->Fill();
        }
    }
    // Write and close the output files
    for (int i = 0; i < 18; ++i) {
        output_files[i]->cd();
        output_trees[i]->Write();
        output_files[i]->Close();
        delete output_files[i];
    }

    input_file->Close();
    delete input_file;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_root_file>" << std::endl;
        return 1;
    }
    const char* input_filename = argv[1];
    process_file(input_filename);

    return 0;
}