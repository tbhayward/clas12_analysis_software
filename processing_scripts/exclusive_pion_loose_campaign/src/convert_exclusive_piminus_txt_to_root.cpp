#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>
#include <TKey.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <stdexcept>
#include <memory>
#include "beam_pol.h"
#include "RunMetadata.h"
#include "LooseRow.h"
using namespace std;
const string schema="runnum evnum helicity beam_energy p1_index p2_index detector1 detector2 sector_e sector1 sector2 fiducial_status e_p e_theta e_phi vz_e p1_p p1_theta p1_phi vz_p1 p2_p p2_theta p2_phi vz_p2 chi2pid1 chi2pid2 beta1 beta2 Q2 W x y nu t tmin tprime Mx2 Mx2_1 Mx2_2 Emiss pTmiss Mh z2 xF2 pT2 phi1 phi2 Delta_phi open_angle_ep1 open_angle_ep2 DepA DepB DepC DepV DepW e_index detector_e status_e status1 status2 e_raw_px e_raw_py e_raw_pz p1_raw_px p1_raw_py p1_raw_pz p2_raw_px p2_raw_py p2_raw_pz chi2pid_e beta_e pass_e_vertex pass_e_sampling pass_e_diagonal pass_e_pcal1 pass_e_pcal2 pass_e_dc pass_p1_momentum pass_p1_theta pass_p1_vertex pass_p1_dc pass_p1_cvt1 pass_p1_cvt2 pass_p1_chi2pid pass_p2_pid pass_p2_vertex pass_p2_dc pass_p2_cvt1 pass_p2_cvt2 pass_e_legacy pass_p1_legacy pass_p2_legacy electron_correction proton_correction pass_dis pass_q2_1p5 kinematics_valid";
vector<string> split(const string& s){istringstream in(s);vector<string> v;string x;while(in>>x)v.push_back(x);return v;}
void validate(const char* path,long long expected){unique_ptr<TFile> f(TFile::Open(path));if(!f||f->IsZombie()||f->TestBit(TFile::kRecovered))throw runtime_error("Invalid ROOT file");auto t=dynamic_cast<TTree*>(f->Get("PhysicsEvents"));if(!t||t->GetEntries()!=expected)throw runtime_error("Tree entry mismatch");auto names=split(schema+" beam_pol target_pol");if(t->GetListOfBranches()->GetEntries()!=int(names.size()))throw runtime_error("Unexpected branch count");for(auto& n:names)if(!t->GetBranch(n.c_str()))throw runtime_error("Missing branch: "+n);}
int main(int argc,char** argv){try{
 if(argc==4&&string(argv[1])=="--validate"){validate(argv[2],stoll(argv[3]));return 0;}
 if(argc!=6)throw runtime_error("Usage: converter input.txt output.root run_info.csv target expected_run; or --validate file.root entries");
 int expectedRun=stoi(argv[5]);string target=argv[4];
 double target_pol=targetPolarization(argv[3],target,expectedRun),beam_pol=getPol(expectedRun);
 string line;ifstream in(argv[1]);if(!in||!getline(in,line)||line!="# "+schema)throw runtime_error("Text schema mismatch");auto names=split(schema);vector<double> d(names.size());vector<int> ints(names.size());set<string> integerNames;for(auto& n:split("runnum evnum helicity p1_index p2_index detector1 detector2 sector_e sector1 sector2 fiducial_status e_index detector_e status_e status1 status2 pass_e_vertex pass_e_sampling pass_e_diagonal pass_e_pcal1 pass_e_pcal2 pass_e_dc pass_p1_momentum pass_p1_theta pass_p1_vertex pass_p1_dc pass_p1_cvt1 pass_p1_cvt2 pass_p1_chi2pid pass_p2_pid pass_p2_vertex pass_p2_dc pass_p2_cvt1 pass_p2_cvt2 pass_e_legacy pass_p1_legacy pass_p2_legacy electron_correction proton_correction pass_dis pass_q2_1p5 kinematics_valid"))integerNames.insert(n);
 TFile out(argv[2],"RECREATE");if(out.IsZombie())throw runtime_error("Cannot create ROOT file");TTree tree("PhysicsEvents","Exclusive e' p pi-; stationary neutron reference");
 for(size_t i=0;i<names.size();i++){string leaf=names[i]+(integerNames.count(names[i])?"/I":"/D");if(integerNames.count(names[i]))tree.Branch(names[i].c_str(),&ints[i],leaf.c_str());else tree.Branch(names[i].c_str(),&d[i],leaf.c_str());}
 tree.Branch("beam_pol",&beam_pol);tree.Branch("target_pol",&target_pol);long long rows=0;
 while(getline(in,line)){if(line.empty())throw runtime_error("Empty data row");parseLooseRow(line,names,d,ints,integerNames,expectedRun);if(tree.Fill()<0)throw runtime_error("ROOT Fill failed");rows++;}
 if(!in.eof())throw runtime_error("Input read failed");if(tree.Write()<=0)throw runtime_error("ROOT Write failed");TNamed convention("kinematic_convention","stationary neutron; t=(q-pi)^2; exact gamma*n->pi-p forward tmin; tprime=t-tmin; radians; opening angles degrees");convention.Write();tree.SetDirectory(nullptr);out.Close();validate(argv[2],rows);cout<<rows<<endl;return 0;
 }catch(const exception& e){cerr<<e.what()<<endl;return 1;}}
